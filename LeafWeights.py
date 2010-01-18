from optparse import OptionParser
import copy
import unittest

import Newick
import FelTree

# A tree like this is used as an example in a manuscript by Eric Stone.
stone_example_tree = '(a:2, (b:1, c:1, d:1, e:1)x:1)y;'

# This tree is from page 125 of Biological Sequence Analysis.
# The Thompson weights for this tree are supposedly
#        w = (20, 20, 32, 47)
# but I have not checked this, so it is left out of the unit tests.
bsa_example_tree = '(((a:2, b:2):3, c:5):3, d:8);'

# Weighting in sequence space: A comparison of methods in terms of generalized sequences.
# Vingron and Sibbald
# Figure 3 and table 3.
# I think the paper gives the wrong weights.
g_vingron_additive_tree = '((((1:7, 2:7):19, (7:16, 8:18):3):12, (5:10, 6:16):28):16, ((3:18, 4:14):15, (9:8, 10:14):24):11);'
g_vingron_ultra_tree = '((((1:6, 2:6):12.9, (7:17, 8:17):1.9):5.3, (5:8.5, 6:8.5):15.7):9.6, ((3:15, 4:15):12.1, (9:11, 10:11):16.1):6.7);'

# Weights for data related by a tree.
# Altschul and Carroll and Lipman, 1989 J. Mol. Biol.
# http://faculty.washington.edu/jht/GS554/Papers/PhylogeneticTrees/AltschulTreeWeights.pdf
# Figure 2 and table 1.
# The rounding of BRU is changed here;
# the authors had two incompatible goals: correct rounding and summing weights to 1.0.
g_acl_tree = '(((((((((BH8:0.7, PV22:0.3, BH10:0.3):0.1, BRU:0.5):0.1, HXB:0.7):2.4, SF2:3.3):0.1, CDC:3.7):0.5, (WMJ1:0.8, WMJ2:0.9, WMJ3:0.9):2.1):0.4, RF:4.3):2.6, ((Z6:2.2, ELI:4.2):2.1, MAL:6.1):1.9):2.7, Z3:9.3);'
g_acl_ordered_names = [
        'BH8', 'PV22', 'BH10', 'BRU', 'HXB',
        'SF2', 'CDC', 'WMJ1', 'WMJ2', 'WMJ3',
        'RF', 'Z6', 'ELI', 'MAL', 'Z3']
g_acl_expected_weights = {
        'BH8' : '0.006',
        'PV22' : '0.013',
        'BH10' : '0.013',
        'BRU' : '0.015',
        'HXB' : '0.017',
        'SF2' : '0.050',
        'CDC' : '0.048',
        'WMJ1' : '0.039',
        'WMJ2' : '0.035',
        'WMJ3' : '0.035',
        'RF' : '0.085',
        'Z6' : '0.129',
        'ELI' : '0.068',
        'MAL' : '0.115',
        'Z3' : '0.333'}



def _add_subtree_resistances(tree):
    """
    This dynamic programming algorithm decorates the directed branches with subtree resistances.
    @param tree: a FelTree object with split branches
    """
    for source, directed_branch in tree.gen_postorder_exits():
        target = directed_branch.get_target()
        if target.is_tip():
            # handle the base case of the recursion
            directed_branch.resistance = 0.0
        else:
            # handle the inductive case of the recursion
            subtree_resistances = []
            for potential_branch in target.gen_directed_branches():
                if potential_branch.get_target() is not source:
                    branch_resistance = potential_branch.get_undirected_branch().get_branch_length()
                    sub_subtree_resistance = potential_branch.resistance
                    subtree_resistances.append(branch_resistance + sub_subtree_resistance)
            directed_branch.resistance = 1.0 / sum(1.0 / r for r in subtree_resistances)

def _add_weights(tree):
    """
    This dynamic programming algorithm decorates the directed branches with weights.
    @param tree: a FelTree object with split branches and subtree resistances and weight allocation for midpoint nodes
    """
    # First initialize all of the weights to zero.
    for source, directed_branch in tree.gen_preorder_exits():
        directed_branch.weight = 0.0
    # Now run the dynamic programming algorithm.
    for source, directed_branch in tree.gen_preorder_exits():
        target = directed_branch.get_target()
        # Initialize the pool of weight to the amount that has been pushed here already.
        weight = directed_branch.weight
        # No weight can be pushed beyond the tips of the tree.
        if target.is_tip():
            target.weight = weight
            continue
        # Do different things depending on whether or not the target is a midpoint.
        if target.weight_allocation is not None:
            # Get the forward and backward branches.
            backward_branch = target.get_directed_branch_to(source)
            forward_branch = [branch for branch in target.gen_directed_branches() if branch is not backward_branch][0]
            # Get the forward and backward resistances.
            backward_resistance = backward_branch.get_undirected_branch().get_branch_length() + backward_branch.resistance
            forward_resistance = forward_branch.get_undirected_branch().get_branch_length() + forward_branch.resistance
            # Get the fraction of weight pushed forward.
            backward_conductance = 1.0 / backward_resistance
            forward_conductance = 1.0 / forward_resistance
            forward_weight_fraction = forward_conductance / (forward_conductance + backward_conductance)
            # Add this new weight to the pool.
            weight += forward_weight_fraction * target.weight_allocation
            # Because we know we are at a midpoint, all of the pool gets pushed straight through.
            forward_branch.weight += weight
        else:
            # Push all existing weight on the directed branch forward through the next node and out the other side.
            target_exits = []
            for potential_branch in target.gen_directed_branches():
                if potential_branch.get_target() is not source:
                    target_exits.append(potential_branch)
            # Get the conductance of each allowed exit from the target.
            exit_conductance_pairs = []
            for branch in target_exits:
                conductance = 1.0 / (branch.get_undirected_branch().get_branch_length() + branch.resistance)
                exit_conductance_pairs.append((branch, conductance))
            # Get the sum of conductances of allowed exits from the target.
            total_conductance = sum(conductance for branch, conductance in exit_conductance_pairs)
            # Push the pool of weight out the exits.
            for branch, conductance in exit_conductance_pairs:
                exit_weight_fraction = conductance / total_conductance
                branch.weight += weight * exit_weight_fraction


def get_stone_weights_fast(tree):
    """
    This is supposed to be an implementation of the weighting of Stone and Sidow.
    It is supposed to have running time on the order of the number of leaves,
    rather than on the order of the square of the number of leaves.
    @param tree: a FelTree object
    @return: a sequence of (name, weight) pairs
    """
    # get the total length of all branches on the tree
    tree_length = tree.get_length()
    # split the branches in half, marking nodes added as midpoints
    old_node_id_set = set(id(node) for node in tree.preorder())
    tree.split_branches()
    for node in tree.preorder():
        if id(node) in old_node_id_set:
            node.weight_allocation = None
        else:
            branch_length = node.get_directed_branch_to_parent().get_undirected_branch().get_branch_length() * 2.0
            node.weight_allocation = branch_length / tree_length
    # calculate subtree resistances using dynamic programming
    _add_subtree_resistances(tree)
    # calculate directed weights using dynamic programming
    _add_weights(tree)
    # the tips are now annotated with weights
    return [(tip.name, tip.weight) for tip in tree.gen_tips()]

def get_thompson_weights(tree):
    """
    Thompson, J. D., Higgins, D. G. and Gibson, T. J. 1994.
    Improved sensitivity of profile searches through the use of sequence weights and gap excision.
    Computer Applications in the Biosciences 10:19-29.
    @param tree: a tree object with branch lengths for all non-root nodes
    @return: a sequence of (name, weight) pairs
    """
    # augment the nodes with subtree resistances calculated postorder
    for node in tree.postorder():
        if node.children:
            node.subtree_resistance = 1.0 / sum(1.0/(child.blen + child.subtree_resistance) for child in node.children)
        else:
            node.subtree_resistance = 0
    # augment the nodes with currents calculated preorder
    for node in tree.preorder():
        if node is tree.root:
            node.current = 1.0
        for child in node.children:
            child.current = (node.current * node.subtree_resistance) / (child.blen + child.subtree_resistance)
    return [(tip.name, tip.current) for tip in tree.gen_tips()]

def get_stone_weights(tree):
    """
    This method was proposed by Stone and Sidow.
    @param tree: a tree object with branch lengths for all non-root nodes
    @return: a sequence of (name, weight) pairs
    """
    # augment each node with an identifier that will survive a deep copy
    for i, node in enumerate(tree.preorder()):
        node.id = i
    # average over all rootings of the tree
    tip_id_to_weight = {}
    for old_target in tree.gen_non_root_nodes():
        # create a new rerooted tree
        clone = copy.deepcopy(tree)
        new_target_list = [node for node in clone.preorder() if node.id == old_target.id]
        assert len(new_target_list) == 1
        target = new_target_list[0]
        new_root = Newick.NewickNode()
        clone.insert_node(new_root, target.parent, target, .5)
        clone.reroot(new_root)
        # find the weights of the rerooted tree using a more traditional method
        # the 'current' attribute added to each tip is its weight
        get_thompson_weights(clone)
        # for each tip add the contribution of this weighting
        for tip in clone.gen_tips():
            weight = tip_id_to_weight.get(tip.id, 0)
            contribution = old_target.blen * tip.current
            tip_id_to_weight[tip.id] = weight + contribution
    # report the final weights
    grand_total_weight = sum(tip_id_to_weight.values())
    return [(tip.name, tip_id_to_weight[tip.id] / grand_total_weight) for tip in tree.gen_tips()]


class TestLeafWeights(unittest.TestCase):

    def do_stone_tree_comparison(self, tree, method, expected_name_weight_pairs):
        """
        Use the example tree by Eric Stone to compare the output of a method to its expected output.
        @param tree: a newick tree of some sort
        @param method: a function that generates (name, weight) pairs
        @param expected_name_weight_pairs: the expected (name, weight) pairs
        """
        names = [tip.name for tip in tree.gen_tips()]
        n = float(len(list(tree.gen_tips())))
        # get the dictionary that maps each name to its calculated weight
        name_to_weight = dict(method(tree))
        # get the dictionary that maps each name to its expected weight
        name_to_expected_weight = dict(expected_name_weight_pairs)
        # for each name compare the calculated weight to the expected weight
        for name in names:
            self.assertAlmostEqual(name_to_weight[name], name_to_expected_weight[name])

    def test_stone_fast(self):
        n = 5.0
        tree = Newick.parse(stone_example_tree, FelTree.NewickTree)
        expected_name_weight_pairs = []
        expected_first_value = (5*n-2) / ((n+2)*(3*n-2))
        expected_non_first_value = (3*n+2) / ((n+2)*(3*n-2))
        expected_name_weight_pairs.append(('a', expected_first_value))
        for name in list('bcde'):
            expected_name_weight_pairs.append((name, expected_non_first_value))
        self.do_stone_tree_comparison(tree, get_stone_weights_fast, expected_name_weight_pairs)

    def test_stone(self):
        n = 5.0
        tree = Newick.parse(stone_example_tree, Newick.NewickTree)
        expected_name_weight_pairs = []
        expected_first_value = (5*n-2) / ((n+2)*(3*n-2))
        expected_non_first_value = (3*n+2) / ((n+2)*(3*n-2))
        expected_name_weight_pairs.append(('a', expected_first_value))
        for name in list('bcde'):
            expected_name_weight_pairs.append((name, expected_non_first_value))
        self.do_stone_tree_comparison(tree, get_stone_weights, expected_name_weight_pairs)

    def test_thompson(self):
        n = 5.0
        tree = Newick.parse(stone_example_tree, Newick.NewickTree)
        expected_name_weight_pairs = []
        expected_first_value = n / (3*n - 2)
        expected_non_first_value = 2 / (3*n - 2)
        expected_name_weight_pairs.append(('a', expected_first_value))
        for name in list('bcde'):
            expected_name_weight_pairs.append((name, expected_non_first_value))
        self.do_stone_tree_comparison(tree, get_thompson_weights, expected_name_weight_pairs)

    def test_stone_fast_equality(self):
        """
        Assert that the faster algorithm gives the same result as the slower algorithm on a nontrivial tree.
        """
        # define the common tree string
        tree_string = bsa_example_tree
        # get (name, weight) pairs using a slow algorithm
        tree = Newick.parse(tree_string, Newick.NewickTree)
        stone_name_weight_pairs = get_stone_weights(tree)
        # get (name, weight) pairs using a faster algorithm
        tree = Newick.parse(tree_string, FelTree.NewickTree)
        fast_stone_name_weight_pairs = get_stone_weights_fast(tree)
        # assert that the algorithms give the same results
        names = [tip.get_name() for tip in tree.gen_tips()]
        da = dict(stone_name_weight_pairs)
        db = dict(fast_stone_name_weight_pairs)
        for name in names:
            self.assertAlmostEqual(da[name], db[name])

    def test_vingron_additive(self):
        tree = Newick.parse(g_vingron_additive_tree, Newick.NewickTree)
        name_weight_pairs = get_thompson_weights(tree)
        int_weight_pairs = [(int(name), weight) for name, weight in name_weight_pairs]
        weights = [weight for name, weight in sorted(int_weight_pairs)]
        #print 'vingron additive:'
        #print weights
        pass

    def test_vingron_ultra(self):
        tree = Newick.parse(g_vingron_ultra_tree, Newick.NewickTree)
        name_weight_pairs = get_thompson_weights(tree)
        int_weight_pairs = [(int(name), weight) for name, weight in name_weight_pairs]
        weights = [weight for name, weight in sorted(int_weight_pairs)]
        #print 'vingron ultra:'
        #print weights
        pass

    def test_acl(self):
        tree = Newick.parse(g_acl_tree, Newick.NewickTree)
        name_weight_pairs = get_thompson_weights(tree)
        name_to_weight = dict(name_weight_pairs)
        self.assertEqual(set(n for n, w in name_weight_pairs), set(g_acl_ordered_names))
        incorrect_names = []
        for name in g_acl_ordered_names:
            s_expected = g_acl_expected_weights[name]
            s_observed = '%.3f' % name_to_weight[name]
            self.assertTrue(s_expected, s_observed)


if __name__ == '__main__':
    unittest.main()


