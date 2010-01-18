"""
Do stuff with Felsenstein's phylogenetic contrasts.
The tree object that is appropriate for this module is a
rooted weighted binary phylogenetic xtree.
Construct the contrast matrix from a rooted binary tree.
Construst the rooted binary tree from the contrast matrix.
"""

import unittest
import math

import numpy as np

import Newick
import FelTree
import NewickIO
import Euclid
import MatrixUtil

# This tree is from figure 25.4 in "Inferring Phylogenies".
g_felsenstein_tree_string = '(((a:0.3, b:0.1):0.25, c:0.65):0.2, (d:0.1, e:0.1):0.7);'

# Columns of this matrix are contrasts from figure 25.4.
g_contrast_matrix = np.array([
    [1, 0.25, 0, 1.0/6],
    [-1, 0.75, 0, 0.5],
    [0, -1, 0, 1.0/3],
    [0, 0, 1, -0.5],
    [0, 0, -1, -0.5]])

# This is the (incorrect) vector of contrast variances from figure 25.4.
g_contrast_variances = (0.4, 0.975, 0.2, 1.11666)


class ContrastError(Exception):
    """
    This is raised when a contrast set is inconsistent with a tree.
    """
    pass


def is_negative(x):
    return x < 0

def is_positive(x):
    return x > 0

def assert_contrast_matrix(C, eps=1e-10):
    """
    @param C: a contrast matrix as a numpy array
    @param eps: numerical error allowed when checking contraints
    """
    if len(C.shape) != 2:
        raise MatrixUtil.MatrixError('the contrast matrix should have two dimensions')
    nrows, ncols = C.shape
    if not (nrows > ncols):
        raise MatrixUtil.MatrixError('the contrast matrix should have more rows than columns')
    for contrast in C.T.tolist():
        if not any(x < 0 for x in contrast):
            raise MatrixUtil.MatrixError('each column of the contrast matrix should have at least one element less than zero')
        if not any(x > 0 for x in contrast):
            raise MatrixUtil.MatrixError('each column of the contrast matrix should have at least one element greater than zero')
        if abs(sum(contrast)) > eps:
            raise MatrixUtil.MatrixError('each column of the contrast matrix should sum to zero')


class ReconstructionInfo:
    """
    This is a node used to help reconstruct the tree from the contrasts.
    These nodes define a scaffold upon which the tree is built.
    """
    
    def __init__(self, node):
        self.node = node
        self.variance = None
        self.child_info_pair = []

    def build_subtree(self, current_contrast, subtree_contrasts, neg_weight, ordered_names):
        """
        This is a recursive function that builds the tree.
        @param current_contrast: the contrast of the current node
        @param subtree_contrasts: a set of contrasts defining the subtree
        @param neg_weight: a way of partitioning the variance or None if root
        @param ordered_names: a list of names conformant with the contrast vector
        """
        # build the negative and positive subtrees
        for i, mycmp in enumerate((is_negative, is_positive)):
            indices = frozenset(i for i, x in enumerate(current_contrast) if mycmp(x))
            child_contrasts = []
            child_subtree_contrasts = []
            for contrast in subtree_contrasts:
                nonzero_indices = frozenset(i for i, x in enumerate(contrast) if x)
                if nonzero_indices == indices:
                    child_contrasts.append(contrast)
                elif nonzero_indices < indices:
                    child_subtree_contrasts.append(contrast)
            if len(child_contrasts) > 1:
                raise ContrastError()
            elif len(child_contrasts) == 0 and len(indices) != 1:
                raise ContrastError()
            # create the child node
            child_node = Newick.NewickNode()
            child_node.parent = self.node
            self.node.children.append(child_node)
            child_info = ReconstructionInfo(child_node)
            if child_contrasts:
                # the child node is an internal node
                child_contrast = child_contrasts[0]
                child_neg_weight = _get_neg_weight(current_contrast, child_contrast)
                child_info.build_subtree(child_contrast, child_subtree_contrasts, child_neg_weight, ordered_names)
            else:
                # the child node is a leaf node
                index, = indices
                child_node.name = ordered_names[index]
                child_info.variance = 0
            self.child_info_pair.append(child_info)
        # define the child branch lengths using the variance partition
        a0 = self.child_info_pair[0].variance
        b0 = self.child_info_pair[1].variance
        total_variance = _get_contrast_variance(current_contrast)
        if neg_weight is None:
            a = (total_variance - a0 - b0) / 2
            b = (total_variance - a0 - b0) / 2
        else:
            a, b = _get_branch_lengths(total_variance, (1 - neg_weight), a0, b0)
        alpha = a + a0
        beta = b + b0
        # set the branch lengths of the children
        self.node.children[0].blen = a
        self.node.children[1].blen = b
        # set the variance of the current node
        self.variance = (alpha * beta) / (alpha + beta)


def _get_contrast_variance(contrast):
    """
    This helper function helps to build a tree from a contrast matrix.
    @param contrast: a contrast vector
    """
    pos_indices = [i for i, x in enumerate(contrast) if x > 0]
    return 1 / sum(contrast[i] for i in pos_indices)**2

def _get_neg_weight(parent_contrast, child_contrast):
    """
    This helper function helps to build a tree from a contrast matrix.
    @param parent_contrast: the parent contrast vector
    @param child_contrast: the child contrast vector
    @return: the proportion of negative weight in the child contrast
    """
    neg_indices = [i for i, x in enumerate(child_contrast) if x < 0]
    pos_indices = [i for i, x in enumerate(child_contrast) if x > 0]
    numerator = sum(parent_contrast[i] for i in neg_indices)
    denominator = sum(parent_contrast[i] for i in neg_indices + pos_indices)
    return numerator / denominator

def _get_branch_lengths(total_variance, a_proportion, partial_a, partial_b):
    """
    This helper function helps to build a tree from a contrast matrix.
    @param total_variance: the sum of the variance of subtrees a and b
    @param a_proportion: the proportion of the variance defined by subtree a
    @param partial_a: the amount of variance defined by the children of subtree a
    @param partial_b: the amount of variance defined by the children of subtree b
    @return: the lengths of branches a and b
    """
    ci = total_variance
    cj = a_proportion
    a0 = partial_a
    b0 = partial_b
    # Solve these two equations for a and b:
    # ci == a + a0 + b + b0
    # cj == (a + a0) / (a + a0 + b + b0)
    alpha = ci * cj
    beta = ci * (1 - cj)
    a = alpha - a0
    b = beta - b0
    return a, b

def contrast_matrix_to_tree(C, ordered_names):
    """
    @param C: contrast matrix as a numpy array
    @param ordered_names: leaf names corresponding to rows of C
    @return: a newick tree object
    """
    contrasts = C.T.tolist()
    # partition the contrasts into the ones with and without entries that are zero
    c_with_zero = [c for c in contrasts if 0 in c]
    c_without_zero = [c for c in contrasts if 0 not in c]
    # exactly one contrast should not have any zero element
    assert len(c_without_zero) == 1
    root_contrast = c_without_zero[0]
    # the variance partition is not defined at the root
    neg_weight = None
    root_node = Newick.NewickNode()
    root_info = ReconstructionInfo(root_node)
    root_info.build_subtree(root_contrast, c_with_zero, neg_weight, ordered_names)
    # get a newick tree from the newick root
    tree = Newick.NewickTree(root_node)
    return tree


class NodeInfo:

    def __init__(self, node):
        """
        @param node: the node of a tree
        """
        self.node = node
        self.resistance = 0
        self.name_weight_pairs = []

    def expand(self, id_to_info):
        """
        Call this only when child nodes have been expanded.
        This requirement could be satisfied by traversing postorder.
        The idea is to avoid deep recursion, using dynamic programming instead.
        @param id_to_info: maps a node id to a NodeInfo object
        """
        self._update_resistance(id_to_info)
        if self.node.has_children():
            for child in self.node.gen_children():
                child_info = id_to_info[id(child)]
                r = child.get_branch_length() + child_info.resistance
                subtree_weight = self.resistance / r
                for name, weight in child_info.name_weight_pairs:
                    self.name_weight_pairs.append((name, subtree_weight * weight))
        else:
            self.name_weight_pairs = [(self.node.get_name(), 1.0)]

    def get_contrast(self, ordered_names, id_to_info):
        """
        Contrasts are zero for names not in the subtree.
        Call this only after the node has been expanded.
        @param ordered_names: the vector of all names
        @param id_to_info: maps a node id to a NodeInfo object
        @return: an ordered list of contrast loadings
        """
        assert self.has_contrast()
        contrast = [0]*len(ordered_names)
        name_to_index = dict((name, i) for i, name in enumerate(ordered_names))
        children = list(self.node.gen_children())
        for child, child_sign in zip(children, (-1.0, 1.0)):
            child_info = id_to_info[id(child)]
            r = child.get_branch_length() + child_info.resistance
            for name, weight in child_info.name_weight_pairs:
                contrast[name_to_index[name]] = weight * child_sign
        return contrast

    def has_contrast(self):
        return self.node.get_child_count() == 2

    def get_variance(self, id_to_info):
        """
        @param id_to_info: maps a node id to a NodeInfo object
        @return: the variance of the contrast that corresponds to this node
        """
        assert self.has_contrast()
        variance = 0
        for child in self.node.gen_children():
            child_info = id_to_info[id(child)]
            variance += child.get_branch_length() + child_info.resistance
        return variance

    def _update_resistance(self, id_to_info):
        """
        @param id_to_info: maps a node id to a NodeInfo object
        """
        if self.node.has_children():
            conductance = 0
            for child in self.node.gen_children():
                child_info = id_to_info[id(child)]
                r = child.get_branch_length() + child_info.resistance
                conductance += 1/r
            self.resistance = 1/conductance
        else:
            self.resistance = 0


def get_contrasts_and_variances(tree, ordered_names):
    """
    @param tree: a newick tree object with branch lengths
    @param ordered_names: a vector of ordered names
    @return: a list of contrasts, and a conformant list of variances
    """
    contrasts = []
    variances = []
    id_to_info = dict((id(node), NodeInfo(node)) for node in tree.preorder())
    for node in tree.postorder():
        info = id_to_info[id(node)]
        info.expand(id_to_info)
        if node.has_children():
            contrast = info.get_contrast(ordered_names, id_to_info)
            variance = info.get_variance(id_to_info)
            contrasts.append(contrast)
            variances.append(variance)
    return contrasts, variances

def get_contrast_matrix(tree, ordered_names):
    """
    @param tree: a newick tree object with branch lengths
    @param ordered_names: a vector of ordered names
    @return: a contrast matrix as a numpy array where columns are contrasts
    """
    contrasts, variances = get_contrasts_and_variances(tree, ordered_names)
    return np.dot(np.array(contrasts).T, np.diag(1/np.sqrt(np.array(variances))))


class TestContrasts(unittest.TestCase):

    def test_felsenstein(self):
        tree = NewickIO.parse(g_felsenstein_tree_string, FelTree.NewickTree)
        ordered_names = ('a', 'b', 'c', 'd', 'e')
        C_expected = np.dot(g_contrast_matrix, np.diag(1/np.sqrt(g_contrast_variances)))
        contrasts, variances = get_contrasts_and_variances(tree, ordered_names)
        C_observed = np.dot(np.array(contrasts).T, np.diag(1/np.sqrt(np.array(variances))))
        """
        print
        print 'felsenstein variances:'
        print g_contrast_variances
        print 'observed variances:'
        print variances
        print
        print 'felsenstein contrast matrix:'
        print C_expected
        print 'observed contrast matrix:'
        print C_observed
        L_expected = np.dot(C_expected, C_expected.T)
        L_observed = np.dot(C_observed, C_observed.T)
        print 'felsenstein L matrix:'
        print L_expected
        print 'observed L matrix:'
        print L_observed
        D = np.array(tree.get_distance_matrix(ordered_names))
        L = Euclid.edm_to_laplacian(D)
        print 'L matrix derived from the D matrix:'
        print L
        """
        pass

    def test_contrast_matrix_to_tree(self):
        original_tree = NewickIO.parse(g_felsenstein_tree_string, FelTree.NewickTree)
        ordered_names = ('a', 'b', 'c', 'd', 'e')
        C = get_contrast_matrix(original_tree, ordered_names)
        assert_contrast_matrix(C)
        reconstructed_tree = contrast_matrix_to_tree(C, ordered_names)
        newick_string = NewickIO.get_newick_string(reconstructed_tree)
        print
        print newick_string
        pass


if __name__ == '__main__':
    unittest.main()
