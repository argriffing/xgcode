"""
Compare a phylogenetic tree to a reference tree.
"""

import unittest

import FelTree
import NewickIO
import Util


def get_nontrivial_split_count(tree):
    """
    @param tree: an unrooted tree object
    @return: the number of nontrivial splits implied by the tree
    """
    return len(get_nontrivial_partitions(tree))

def get_weighted_split_count(tree):
    """
    A balanced nontrivial split adds more weight than a less balanced nontrivial split.
    @param tree: an unrooted tree object
    @return: the weighted number of nontrivial splits implied by the tree
    """
    parts = get_nontrivial_partitions(tree)
    total = 0
    for a, b in parts:
        n = len(a) + len(b)
        k = min(len(a), len(b))
        total += Util.choose(n, k)
    return total

def _get_branch_id_to_node_id_set(tree):
    """
    I want this function so I can get full id splits of the tree.
    By full I mean including internal nodes.
    @param tree: a tree object
    @return: a map from the id fo a directed branch to a set of node ids
    """
    d = {}
    for source, dbranch in tree.gen_postorder_exits():
        id_set = set()
        target = dbranch.get_target()
        id_set.add(id(target))
        for next_dbranch in target.gen_exits(source):
            id_set.update(d[id(next_dbranch)])
        d[id(dbranch)] = id_set
    return d

def _get_branch_id_to_leaf_name_set(tree):
    """
    @param tree: a tree object
    @return: a map from the id of a directed branch to a set of leaf names
    """
    d = {}
    for source, dbranch in tree.gen_postorder_exits():
        leaf_name_set = set()
        target = dbranch.get_target()
        if target.is_tip():
            leaf_name_set.add(target.get_name())
        else:
            for next_dbranch in target.gen_exits(source):
                leaf_name_set.update(d[id(next_dbranch)])
        d[id(dbranch)] = leaf_name_set
    return d

def get_partitions(tree):
    """
    Get all of the partitions implied by a tree.
    Each positive branch implies a partition.
    Each partition is a frozenset of two frozensets of leaf names.
    The return value is the set of these partitions.
    Note that the word 'partition' is a python keyword,
    so the word 'part' will be used here instead.
    @param tree: a tree object
    @return: the set of partitions implied by the tree.
    """
    # map a directed branch id to the set of leaf names in its subtree
    d = _get_branch_id_to_leaf_name_set(tree)
    # for each branch in the tree get the frozenset of leaf names on each end of the branch
    parts = set()
    for node in tree.gen_non_root_nodes():
        parent = node.get_parent()
        directed_branches = (node.get_directed_branch_to(parent), parent.get_directed_branch_to(node))
        branch_length = directed_branches[0].get_branch_length()
        for dbranch in directed_branches:
            assert dbranch.get_branch_length() == branch_length
        if branch_length > 0:
            leaf_sets = [d[id(dbranch)] for dbranch in directed_branches]
            part = frozenset(frozenset(leaf_set) for leaf_set in leaf_sets)
            parts.add(part)
    # return the set of partitions
    return parts

def get_partitions_and_branch_lengths(tree):
    """
    Each partition is a frozenset of two frozensets of leaf names.
    @param tree: a tree object
    @return: a set of (partition, branch length) pairs
    """
    # map a directed branch id to the set of leaf names in its subtree
    d = _get_branch_id_to_leaf_name_set(tree)
    # for each branch in the tree get the frozenset of leaf names on each end of the branch
    ret = set()
    for node in tree.gen_non_root_nodes():
        parent = node.get_parent()
        directed_branches = (node.get_directed_branch_to(parent), parent.get_directed_branch_to(node))
        branch_length = directed_branches[0].get_branch_length()
        for dbranch in directed_branches:
            assert dbranch.get_branch_length() == branch_length
        if branch_length > 0:
            leaf_sets = [d[id(dbranch)] for dbranch in directed_branches]
            part = frozenset(frozenset(leaf_set) for leaf_set in leaf_sets)
            ret.add((part, branch_length))
    # return the set of (partition, branch length) pairs
    return ret

def get_nontrivial_partitions(tree):
    """
    Get all of the nontrivial partitions implied by a tree.
    Note that the word 'partition' is a python keyword,
    so the word 'part' will be used here instead.
    @param tree: a tree object
    @return: the set of nontrivial partitions implied by the tree.
    """
    nontrivial_partitions = set()
    for part in get_partitions(tree):
        if min(len(leaf_set) for leaf_set in part) > 1:
            nontrivial_partitions.add(part)
    return nontrivial_partitions

def get_split_distance(observed_tree, expected_tree):
    """
    @param observed_tree: a tree object
    @param expected_tree: a tree object
    @return: the number of nontrivial splits in the expected tree that are not in the observed tree
    """
    expected_partitions = get_nontrivial_partitions(expected_tree)
    observed_partitions = get_nontrivial_partitions(observed_tree)
    missing_partitions = expected_partitions - observed_partitions
    return len(missing_partitions)

def get_weighted_split_distance(observed_tree, expected_tree):
    """
    @param observed_tree: a tree object
    @param expected_tree: a tree object
    @return: the cost of splits in the expected tree that are not present in the observed tree
    """
    missing_partitions = get_nontrivial_partitions(expected_tree) - get_nontrivial_partitions(observed_tree)
    # get the sum of the count pair weights
    total = 0
    for a, b in missing_partitions:
        n = len(a) + len(b)
        k = min(len(a), len(b))
        total += Util.choose(n, k)
    return total

def _make_partition(first_group, second_group):
    """
    This is a helper function that returns a partition.
    @param first_group: a collection of hashable items
    @param second_group: another collection of hashable items
    @return: a frozenset of two frozensets
    """
    return frozenset((frozenset(first_group), frozenset(second_group)))


class TestTreeComparison(unittest.TestCase):

    def test_get_nontrivial_split_count(self):
        """
        Test the function that gets the number of nontrivial splits
        """
        # define some trees
        tree_string_a = '((A:1, B:1):1, (C:1, D:1):1, (E:1, F:1):1);'
        tree_string_b = '(((A:1, B:1):1, C:1):1, D:1, (E:1, F:1):1);'
        tree_string_c = '((A:1, B:1, C:1):1, D:1, (E:1, F:1):1);'
        tree_string_d = '(((A:1, B:1):1, C:1):1, (D:1, (E:1, F:1):1):1);'
        tree_a = NewickIO.parse(tree_string_a, FelTree.NewickTree)
        tree_b = NewickIO.parse(tree_string_b, FelTree.NewickTree)
        tree_c = NewickIO.parse(tree_string_c, FelTree.NewickTree)
        tree_d = NewickIO.parse(tree_string_d, FelTree.NewickTree)
        # assert that the correct split count is recovered
        self.assertEqual(get_nontrivial_split_count(tree_a), 3)
        self.assertEqual(get_nontrivial_split_count(tree_b), 3)
        self.assertEqual(get_nontrivial_split_count(tree_c), 2)
        self.assertEqual(get_nontrivial_split_count(tree_d), 3)

    def test_get_weighted_split_count(self):
        """
        Test the function that gets the weighted number of nontrivial splits
        """
        # define some trees
        tree_string_a = '((A:1, B:1):1, (C:1, D:1):1, (E:1, F:1):1);'
        tree_string_b = '(((A:1, B:1):1, C:1):1, D:1, (E:1, F:1):1);'
        tree_string_c = '(((A:1, B:1):1, C:1):1, (D:1, (E:1, F:1):1):1);'
        tree_a = NewickIO.parse(tree_string_a, FelTree.NewickTree)
        tree_b = NewickIO.parse(tree_string_b, FelTree.NewickTree)
        tree_c = NewickIO.parse(tree_string_c, FelTree.NewickTree)
        # the weighted split counts are different,
        # even though both trees have internal nodes of order 3 and have the same number of leaves
        self.assertEqual(get_weighted_split_count(tree_a), 45)
        self.assertEqual(get_weighted_split_count(tree_b), 50)
        self.assertEqual(get_weighted_split_count(tree_c), 50)

    def test_get_partitions(self):
        """
        Test the function that gets the set of partitions implied by a tree.
        """
        # get the observed partitions
        tree_string = '((A:1, B:1):1, (C:1, (D:1, E:1):1):1);'
        tree = NewickIO.parse(tree_string, FelTree.NewickTree)
        observed_partitions = get_partitions(tree)
        # get the expected partitions
        arr = (
                ('A', 'BCDE'),
                ('B', 'ACDE'),
                ('C', 'ABDE'),
                ('D', 'ABCE'),
                ('E', 'ABCD'),
                ('AB', 'CDE'),
                ('ABC', 'DE')
                )
        expected_partitions = set()
        for a, b in arr:
            part = frozenset([frozenset(a), frozenset(b)])
            expected_partitions.add(part)
        # assert that the observed partitions equal the expected partitions
        self.assertEqual(observed_partitions, expected_partitions)

    def test_get_nontrivial_partitions(self):
        """
        Test the function that gets the set of nontrivial partitions implied by a tree.
        """
        # get the observed nontrivial partitions
        tree_string = '((A:1, B:1):1, (C:1, (D:1, E:1):1):1);'
        tree = NewickIO.parse(tree_string, FelTree.NewickTree)
        observed_nontrivial_partitions = get_nontrivial_partitions(tree)
        # get the expected nontrivial partitions
        arr = (
                ('AB', 'CDE'),
                ('ABC', 'DE')
                )
        expected_nontrivial_partitions = set()
        for a, b in arr:
            part = frozenset([frozenset(a), frozenset(b)])
            expected_nontrivial_partitions.add(part)
        # assert that the observed nontrivial partitions equal the expected nontrivial partitions
        self.assertEqual(observed_nontrivial_partitions, expected_nontrivial_partitions)
        # each nontrivial partition corresponds to a split
        nontrivial_split_count = get_nontrivial_split_count(tree)
        self.assertEqual(nontrivial_split_count, len(expected_nontrivial_partitions))

    def test_get_nontrivial_partitions_b(self):
        """
        Test nontrivial partitions of trees with branch lengths.
        """
        tree_string_a = '(c:1, (a:1, e:1):1, ((d:1, f:1):1, b:1):1);'
        tree_string_b = '((f:1, d:1):1, ((b:1, (c:1, a:1):1):1, e:1):1);'
        tree_a = NewickIO.parse(tree_string_a, FelTree.NewickTree)
        tree_b = NewickIO.parse(tree_string_b, FelTree.NewickTree)
        partitions_a = get_nontrivial_partitions(tree_a)
        partitions_b = get_nontrivial_partitions(tree_b)
        expected_a = set((
                _make_partition('ae', 'cdfb'),
                _make_partition('aec', 'dfb'),
                _make_partition('aecb', 'df')))
        expected_b = set((
                _make_partition('ac', 'bedf'),
                _make_partition('acb', 'edf'),
                _make_partition('acbe', 'df')))
        self.assertEqual(partitions_a, expected_a)
        self.assertEqual(partitions_b, expected_b)

    def test_get_nontrivial_partitions_c(self):
        """
        Test nontrivial partitions of trees without branch lengths.
        """
        """
        tree_string_a = '(c, (a, e), ((d, f), b));'
        tree_string_b = '((f, d), ((b, (c, a)), e));'
        tree_a = NewickIO.parse(tree_string_a, FelTree.NewickTree)
        tree_b = NewickIO.parse(tree_string_b, FelTree.NewickTree)
        partitions_a = get_nontrivial_partitions(tree_a)
        partitions_b = get_nontrivial_partitions(tree_b)
        expected_a = set((
                _make_partition('ae', 'cdfb'),
                _make_partition('aec', 'dfb'),
                _make_partition('aecb', 'df')))
        expected_b = set((
                _make_partition('ac', 'bedf'),
                _make_partition('acb', 'edf'),
                _make_partition('acbe', 'df')))
        self.assertEqual(partitions_a, expected_a)
        self.assertEqual(partitions_b, expected_b)
        """
        pass

    def test_get_split_distance(self):
        """
        Test the function that gets the number of missing nontrivial partitions.
        """
        # define some trees
        tree_string_a = '((A:1, B:1):1, C:1, (D:1, E:1):1);'
        tree_string_b = '((A:1, B:1):1, D:1, (C:1, E:1):1);'
        tree_string_c = '((A:1, D:1):1, C:1, (B:1, E:1):1);'
        tree_string_d = '((A:1, D:1):1, (C:1, B:1, E:1):1);'
        tree_a = NewickIO.parse(tree_string_a, FelTree.NewickTree)
        tree_b = NewickIO.parse(tree_string_b, FelTree.NewickTree)
        tree_c = NewickIO.parse(tree_string_c, FelTree.NewickTree)
        tree_d = NewickIO.parse(tree_string_d, FelTree.NewickTree)
        # the distance from a tree to itself should be zero
        self.assertEqual(get_split_distance(tree_a, tree_a), 0)
        self.assertEqual(get_split_distance(tree_b, tree_b), 0)
        self.assertEqual(get_split_distance(tree_c, tree_c), 0)
        self.assertEqual(get_split_distance(tree_d, tree_d), 0)
        # some of the distances are symmetric
        self.assertEqual(get_split_distance(tree_a, tree_b), 1)
        self.assertEqual(get_split_distance(tree_b, tree_a), 1)
        self.assertEqual(get_split_distance(tree_b, tree_c), 2)
        self.assertEqual(get_split_distance(tree_c, tree_b), 2)
        self.assertEqual(get_split_distance(tree_a, tree_c), 2)
        self.assertEqual(get_split_distance(tree_c, tree_a), 2)
        # it is possible for the distance to be asymmetric if internal nodes are not order 3
        self.assertEqual(get_split_distance(tree_a, tree_d), 1)
        self.assertEqual(get_split_distance(tree_d, tree_a), 2)

    def test_get_weighted_split_distance(self):
        """
        Test the function that gets the number of missing nontrivial partitions.
        """
        # define some trees
        tree_string_a = '((A:1, B:1):1, (C:1, D:1):1, (E:1, F:1):1);'
        tree_string_b = '(((A:1, B:1):1, C:1):1, D:1, (E:1, F:1):1);'
        tree_a = NewickIO.parse(tree_string_a, FelTree.NewickTree)
        tree_b = NewickIO.parse(tree_string_b, FelTree.NewickTree)
        # the distance from a tree to itself should be zero
        self.assertEqual(get_weighted_split_distance(tree_a, tree_a), 0)
        self.assertEqual(get_weighted_split_distance(tree_b, tree_b), 0)
        # the distance is not necessarily symmetric
        self.assertEqual(get_weighted_split_distance(tree_a, tree_b), 20)
        self.assertEqual(get_weighted_split_distance(tree_b, tree_a), 15)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestTreeComparison)
    unittest.TextTestRunner(verbosity=2).run(suite)

