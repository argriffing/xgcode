"""Construct an example tree with pedagogically useful properties.
"""

from StringIO import StringIO
import time
import random
import optparse

import numpy as np

from SnippetUtil import HandlingError
import MatrixUtil
import SchurAlgebra
import Euclid
import TreeSampler
import BranchLenSampler
import BuildTreeTopology
import Dendrogram
import Xtree
from Form import RadioItem
from Form import CheckItem
import Form
import FormOut

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.RadioGroup('leaf_options', 'number of leaves', [
                RadioItem('six_leaves', '6', True),
                RadioItem('seven_leaves', '7')]),
            Form.CheckGroup('requirements', 'requirements', [
                CheckItem('invalid_dendrogram',
                    'the topology of the naive dendrogram must be incorrect',
                    True),
                CheckItem('informative_children',
                    'splits of the child trees must be informative',
                    True),
                CheckItem('force_difference',
                    'full graph split must differ from Schur graph split',
                    False),
                CheckItem('informative_full_split',
                    'full graph split must be informative', False)]),
            Form.CheckGroup('sampling_options', 'branch length sampling', [
                CheckItem('allow_integers',
                    'allow single digit integer branch lengths', True),
                CheckItem('allow_reciprocals',
                    'allow reciprocal single digit integer branch lengths',
                    False)])]
    return form_objects

def get_form_out():
    return FormOut.Report()


class TreeSearch:
    """
    This is a virtual base class.
    """

    def __init__(self):
        # boolean requirements defined by the user
        self.informative_children = None
        self.force_difference = None
        self.informative_full_split = None
        self.invalid_dendrogram = None
        # search options defined by the subclass
        self.tree = None
        self.desired_primary_split = None
        self.id_to_index = None
        # initialize the counts that are tracked for bookkeeping
        self.aug_split_collision_count = 0
        self.aug_split_degenerate_count = 0
        self.error_primary_split_count = 0
        self.invalid_primary_split_count = 0
        self.degenerate_primary_split_count = 0
        self.undesired_primary_split_count = 0
        self.desired_primary_split_count = 0
        self.uninformative_child_count = 0
        self.informative_child_count = 0
        self.valid_dendrogram_count = 0
        self.success_count = 0

    def is_initialized(self):
        required_data = [
                self.informative_children,
                self.force_difference,
                self.informative_full_split,
                self.invalid_dendrogram,
                self.tree,
                self.desired_primary_split,
                self.id_to_index]
        return not (None in required_data)

    def get_result_text(self):
        """
        @return: a multi-line string of text
        """
        out = StringIO()
        if self.force_difference or self.informative_full_split:
            print >> out, 'full graph split stats:'
            print >> out, self.aug_split_collision_count, 'full graph splits collided with the desired primary split'
            print >> out, self.aug_split_degenerate_count, 'full graph splits were degenerate'
            print >> out
        print >> out, 'primary split stats:'
        print >> out, self.error_primary_split_count, 'errors in finding the primary split (should be 0)'
        print >> out, self.invalid_primary_split_count, 'invalid primary splits (should be 0)'
        print >> out, self.degenerate_primary_split_count, 'degenerate primary splits'
        print >> out, self.undesired_primary_split_count, 'primary splits were not the target split'
        print >> out, self.desired_primary_split_count, 'primary splits were the target split'
        print >> out
        if self.informative_children:
            print >> out, 'secondary split stats:'
            print >> out, self.uninformative_child_count, 'samples had at least one uninformative child tree'
            print >> out, self.informative_child_count, 'samples had two informative child trees'
            print >> out
        if self.invalid_dendrogram:
            print >> out, 'naive dendrogram stats:'
            print >> out, self.valid_dendrogram_count, 'naive dendrograms were valid'
            print >> out
        return out.getvalue().strip()

    def do_search(self, nseconds, sampling_function):
        """
        @param nseconds: the allowed time for the search or None to search until interrupted
        @param sampling_function: a function that samples a branch length
        @return: True if a tree was found that met the criteria
        """
        if not self.is_initialized():
            raise RuntimeError('the search was not sufficiently initialized')
        true_splits = self.tree.get_nontrivial_splits()
        start_time = time.time()
        while True:
            elapsed_time = time.time() - start_time
            if nseconds and elapsed_time > nseconds:
                return False
            # assign new sampled branch lengths
            for branch in self.tree.get_branches():
                branch.length = sampling_function()
            # get the distance matrix so we can use a library function to get the split
            D = np.array(self.tree.get_distance_matrix())
            ntips = len(D)
            # get the Laplacian matrix of the full tree and the corresponding Fiedler split of the leaves
            if self.force_difference or self.informative_full_split:
                A_aug = np.array(self.tree.get_weighted_adjacency_matrix(self.id_to_index))
                L_aug = Euclid.adjacency_to_laplacian(A_aug)
                v_aug = BuildTreeTopology.laplacian_to_fiedler(L_aug)
                left_aug, right_aug = BuildTreeTopology.eigenvector_to_split(v_aug)
                left = [x for x in left_aug if x in range(ntips)]
                right = [x for x in right_aug if x in range(ntips)]
                leaf_eigensplit_aug = BuildTreeTopology.make_split(left, right)
                if self.force_difference:
                    if leaf_eigensplit_aug == self.desired_primary_split:
                        self.aug_split_collision_count += 1
                        continue
                if self.informative_full_split:
                    if min(len(s) for s in leaf_eigensplit_aug) < 2:
                        self.aug_split_degenerate_count += 1
                        continue
            # get the eigensplit
            try:
                eigensplit = BuildTreeTopology.split_using_eigenvector(D)
            except BuildTreeTopology.DegenerateSplitException, e:
                self.degenerate_primary_split_count += 1
                continue
            except BuildTreeTopology.InvalidSpectralSplitException, e:
                self.error_primary_split_count += 1
                continue
            if eigensplit not in true_splits:
                raise RuntimeError('INVALID SPLIT:' + tree.get_newick_string())
            if eigensplit != self.desired_primary_split:
                self.undesired_primary_split_count += 1
                continue
            self.desired_primary_split_count += 1
            # check the splits of the two child trees
            degenerate_subsplit_count = 0
            L = Euclid.edm_to_laplacian(D)
            for side in eigensplit:
                L_child = SchurAlgebra.mmerge(L, side)
                v = BuildTreeTopology.laplacian_to_fiedler(L_child)
                child_eigensplit = BuildTreeTopology.eigenvector_to_split(v)
                if min(len(s) for s in child_eigensplit) < 2:
                    degenerate_subsplit_count += 1
            if degenerate_subsplit_count:
                self.uninformative_child_count += 1
            else:
                self.informative_child_count += 1
            if self.informative_children:
                if degenerate_subsplit_count:
                    continue
            # check the dendrogram
            if self.invalid_dendrogram:
                labels = range(len(D))
                hierarchy = Dendrogram.get_hierarchy(D, Dendrogram.spectral_split, labels)
                dendrogram_splits = set(Dendrogram.hierarchy_to_nontrivial_splits(hierarchy))
                if dendrogram_splits == true_splits:
                    self.valid_dendrogram_count += 1
                    continue
            # the tree has met all of the requirements
            return True


class SevenLeafSearch(TreeSearch):

    def __init__(self):
        """
        Define the topology of a tree for which branch lengths will be sought.
        """
        TreeSearch.__init__(self)
        # create the fixed tree topology
        topo = [[0, 1], [2, 3], [[4, 5], 6]]
        self.tree = Xtree.list_to_uniformly_weighted_tree(topo)
        # define the expected primary split
        self.desired_primary_split = frozenset([frozenset([0, 1, 2, 3]), frozenset([4, 5, 6])])
        # create the mapping from node id to node index
        self.id_to_index = dict((id(node), node.label) for node in self.tree.get_labeled_vertices())
        # define the internal nodes in the left hand subtree
        self.id_to_index[id(self.tree.children[0])] = 7
        self.id_to_index[id(self.tree.children[1])] = 8
        self.id_to_index[id(self.tree)] = 9
        # define the internal nodes in the right hand subtree
        self.id_to_index[id(self.tree.children[2].children[0])] = 10
        self.id_to_index[id(self.tree.children[2])] = 11


class SixLeafSearch(TreeSearch):

    def __init__(self):
        """
        Define the topology of a tree for which branch lengths will be sought.
        """
        TreeSearch.__init__(self)
        # create the fixed tree topology
        topo = [0, [1, 2], [[3, 4], 5]]
        self.tree = Xtree.list_to_uniformly_weighted_tree(topo)
        # define the expected primary split
        self.desired_primary_split = frozenset([frozenset([0, 1, 2]), frozenset([3, 4, 5])])
        # create the mapping from node id to node index
        self.id_to_index = dict((id(node), node.label) for node in self.tree.get_labeled_vertices())
        # define the internal nodes in the left hand subtree
        self.id_to_index[id(self.tree.children[1])] = 6
        self.id_to_index[id(self.tree)] = 7
        # define the internal nodes in the right hand subtree
        self.id_to_index[id(self.tree.children[2].children[0])] = 8
        self.id_to_index[id(self.tree.children[2])] = 9


def do_tree_search(tree_search, nseconds, sampling_function):
    """
    @param tree_search: a tree searching object
    @param nseconds: the allowed time for the search or None to search until interrupted
    @param sampling_function: a function that samples a branch length
    """
    out = StringIO()
    try:
        success = tree_search.do_search(nseconds, sampling_function)
    except KeyboardInterrupt, e:
        success = False
    if success:
        print >> out, 'Found a tree that satisfies the criteria!'
        print >> out, tree_search.tree.get_newick_string()
    else:
        print >> out, 'No tree satisfied the criteria.'
    print >> out
    print >> out, tree_search.get_result_text()
    return out.getvalue().strip()

def get_response_content(fs):
    nseconds = 2
    sampling_function = BranchLenSampler.ShortAscii(
            fs.allow_integers, fs.allow_reciprocals)
    if fs.six_leaves:
        tree_search = SixLeafSearch()
    elif fs.seven_leaves:
        tree_search = SevenLeafSearch()
    tree_search.informative_children = fs.informative_children
    tree_search.force_difference = fs.force_difference
    tree_search.informative_full_split = fs.informative_full_split
    tree_search.invalid_dendrogram = fs.invalid_dendrogram
    return do_tree_search(tree_search, nseconds, sampling_function) + '\n'

def main(options):
    #FIXME use argparse and a nonnegative type instead of the assertion
    assert 0 <= options.nseconds
    sampling_function = BranchLenSampler.ShortAscii(
            options.allow_integers, options.allow_reciprocals)
    tree_search = SixLeafSearch()
    tree_search.informative_children = options.informative_children
    tree_search.force_difference = options.force_difference
    tree_search.informative_full_split = options.informative_full_split
    tree_search.invalid_dendrogram = options.invalid_dendrogram
    print do_tree_search(tree_search, options.nseconds, sampling_function)

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('--nseconds', dest='nseconds', type='int', default=0, help='seconds to run or 0 to run until ctrl-c')
    parser.add_option('--no-reciprocals', action='store_false', dest='allow_reciprocals', default=True, help='disallow reciprocals of single digit branch lengths')
    parser.add_option('--no-integers', action='store_false', dest='allow_integers', default=True, help='disallow single digit branch lengths')
    parser.add_option('--informative-children', action='store_true', dest='informative_children', default=False, help='require informative secondary splits')
    parser.add_option('--force-difference', action='store_true', dest='force_difference', default=False, help='require the Fiedler split to differ from the Schur split')
    parser.add_option('--informative-full-split', action='store_true', dest='informative_full_split', default=False, help='require the Fiedler split to be informative')
    parser.add_option('--invalid-dendrogram', action='store_true', dest='invalid-dendrogram', default=False, help='require the naive dendrogram to be invalid')
    options, args = parser.parse_args()
    main(options)

