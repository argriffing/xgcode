"""
Get all nontrivial splits that would be implied by a tree built from a given distance matrix.
This allows methods other than neighbor joining to be tested.
The idea is to provide two independent components;
first, a function that splits a distance matrix,
second, a function that returns a new distance matrix
given a set of indices to be used as an outgroup.
In neighbor joining this updating function does not affect
pairwise distances between taxa not in the outgroup,
but we want to test new updating methods that could change these distances.
"""

import unittest
import itertools

import numpy as np

import NewickIO
import FelTree
import MatrixUtil
import SchurAlgebra
import Euclid
import Xtree
import iterutils

def eigenvector_to_split(v, epsilon=1e-14):
    """
    This is a helper function.
    @param v: the signs of the loadings of this eigenvector define the split
    @param epsilon: negligible eigenvector loadings will be treated as zero
    @return: a split
    """
    vprime = [0.0 if abs(x) < epsilon else x for x in v]
    positive_indices = [i for i, x in enumerate(vprime) if x > 0]
    nonpositive_indices = [i for i, x in enumerate(vprime) if x <= 0]
    return make_split(positive_indices, nonpositive_indices)

def make_split(a, b):
    """
    This is a helper function.
    @param a: a sequence of hashable values
    @param b: a sequence of hashable values
    @return: a split
    """
    return frozenset((frozenset(a), frozenset(b)))

def get_quartet_score(quartet, D):
    """
    I think this is an estimate of the the internal branch length.
    This is a helper function for evaluating distance matrix criteria.
    @param quartet: a quartet of indices in D
    @param D: a distance matrix
    return: a floating point number
    """
    (i, j), (k, l) = quartet
    return (D[i,k]+D[i,l]+D[j,k]+D[j,l])/2 - D[i,j] - D[k,l]

def leaf_is_interior(x, quartet, quartets):
    """
    See if a leaf is interior to a quartet in a tree.
    This is a helper function for evaluating distance matrix criteria.
    @param x: the integer label of a leaf of an xtree
    @param quartet: a quartet in the tree
    @param quartets: all quartets in the tree
    @return: true if the leaf is interior to the quartet in the tree
    """
    (i, j), (k, l) = quartet
    invalid_quartets = set((
        make_split((i, k), (x, l)),
        make_split((i, k), (x, j)),
        make_split((i, x), (j, l)),
        make_split((k, x), (j, l))))
    # If there is any intersection between the invalid quartets
    # and the true quartets then the leaf is not interior.
    return not (invalid_quartets & quartets)

def is_quartet_additive(tree, D):
    """
    This criterion is mentioned in "Why neighbor joining works".
    @param tree: the root of a weighted phylogenetic xtree
    @param D: a distance matrix conformant to the tree labels
    @return: True if D is quartet additive with respect to the tree
    """
    all_labels = set(tree.get_labels())
    quartets = tree.get_quartets()
    for quartet in quartets:
        (i, j), (k, l) = quartet
        left, right = quartet
        other_labels = all_labels - (left | right)
        for x in other_labels:
            for y in other_labels:
                if x != y:
                    if leaf_is_interior(x, quartet, quartets):
                        if not leaf_is_interior(y, quartet, quartets):
                            qa = make_split((i, j), (x, y))
                            qb = make_split((k, l), (x, y))
                            sa = get_quartet_score(qa, D)
                            sb = get_quartet_score(qb, D)
                            if qa not in quartets:
                                if sb <= sa:
                                    return False
                            if qb not in quartets:
                                if sa <= sb:
                                    return False
    return True

def is_quartet_consistent(tree, D):
    """
    This criterion is mentioned in "Why neighbor joining works".
    @param tree: the root of a weighted phylogenetic xtree
    @param D: a distance matrix conformant to the tree labels
    @return: True if D is quartet consistent with respect to the tree
    """
    for quartet in tree.get_quartets():
        score = get_quartet_score(quartet, D)
        (i,j),(k,l) = quartet
        if score <= get_quartet_score(((i,k),(j,l)), D):
            return False
        if score <= get_quartet_score(((i,l),(j,k)), D):
            return False
    return True

def is_atteson(tree, D):
    """
    The Atteson condition is that each distance is wrong by no more than d/2.
    The value d is the minimum branch length of the tree.
    @param tree: the root of a weighted binary phylogenetic xtree
    @param D: a distance matrix conformant to the tree labels
    @return: True if the Atteson condition holds
    """
    min_branch_length = min(branch.length for branch in tree.get_branches())
    n = len(D)
    D_true = tree.get_distance_matrix()
    for i in range(n):
        for j in range(n):
            abs_difference = abs(D[i][j] - D_true[i][j])
            if abs_difference > min_branch_length / 2.0:
                return False
    return True


class InvalidSpectralSplitException(Exception):
    """
    This exception is raised when a spectral split fails horribly.
    In particular, it is raised when an eigenvector based splitting method
    gets an eigenvector that does not define any split, even a degenerate split.
    """
    def __init__(self, D):
        Exception.__init__(self)
        self.D = D


class DegenerateSplitException(Exception):
    """
    This exception is raised when one taxon is split from the rest.
    It is caught within this module.
    """
    def __init__(self, index):
        Exception.__init__(self)
        self.index = index


def full_split(selected_items, all_items):
    """
    This is a helper function which creates a split from a selection.
    @param selected_items: items on one side of the split
    @param all_items: all of the items to be split
    @return: a frozenset of the two implied disjoint frozensets
    """
    total = set(all_items)
    selection = frozenset(selected_items)
    complement = frozenset(total - selection)
    return frozenset((selection, complement))

def index_split_to_label_split(index_split, label_sets):
    """
    This is a helper function which creates a label split from an index split.
    @param index_split: a split of indices of the label sets
    @param label_sets: a set of labels for each index
    @return: a label split defined by the index split of the label sets
    """
    label_split = set()
    for index_selection in index_split:
        labels = set()
        for i in index_selection:
            labels.update(label_sets[i])
        label_split.add(frozenset(labels))
    return frozenset(label_split)

def get_splits(initial_distance_matrix, split_function, update_function, on_label_split=None):
    """
    This is the most external of the functions in this module.
    Get the set of splits implied by the tree that would be reconstructed.
    @param initial_distance_matrix: a distance matrix
    @param split_function: takes a distance matrix and returns an index split
    @param update_function: takes a distance matrix and an index subset and returns a distance matrix
    @param on_label_split: notifies the caller of the label split induced by an index split
    @return: a set of splits
    """
    n = len(initial_distance_matrix)
    # keep a stack of (label_set_per_vertex, distance_matrix) pairs
    initial_state = ([set([i]) for i in range(n)], initial_distance_matrix)
    stack = [initial_state]
    # process the stack in a depth first manner, building the split set
    label_split_set = set()
    while stack:
        label_sets, D = stack.pop()
        # if the matrix is small then we are done
        if len(D) < 4:
            continue
        # split the indices using the specified function
        try:
            index_split = split_function(D)
            # convert the index split to a label split
            label_split = index_split_to_label_split(index_split, label_sets)
            # notify the caller if a callback is requested
            if on_label_split:
                on_label_split(label_split)
            # add the split to the master set of label splits
            label_split_set.add(label_split)
            # for large matrices create the new label sets and the new conformant distance matrices
            a, b = index_split
            for index_selection, index_complement in ((a,b),(b,a)):
                if len(index_complement) > 2:
                    next_label_sets = SchurAlgebra.vmerge(label_sets, index_selection)
                    next_D = update_function(D, index_selection)
                    next_state = (next_label_sets, next_D)
                    stack.append(next_state)
        except DegenerateSplitException, e:
            # we cannot recover from a degenerate split unless there are more than four indices
            if len(D) <= 4:
                continue
            # with more than four indices we can fall back to partial splits
            index_set = set([e.index])
            # get the next label sets
            next_label_sets = SchurAlgebra.vdelete(label_sets, index_set)
            # get the next conformant distance matrix by schur complementing out the offending index
            L = Euclid.edm_to_laplacian(D)
            L_small = SchurAlgebra.mschur(L, index_set)
            next_D = Euclid.laplacian_to_edm(L_small)
            next_state = (next_label_sets, next_D)
            stack.append(next_state)
    return label_split_set

def split_nj(D):
    """
    Split a distance matrix according to the neighbor joining criterion.
    @param D: a distance matrix
    @return: a set of two index sets defining a split of the indices
    """
    n = len(D)
    # there is no reason to split a 3x3 distance matrix
    assert n > 3, n
    Q = Euclid.edm_to_q(D)
    pairs = [(i, j) for i in range(n) for j in range(n) if i < j]
    best_value, best_pair = min((Q[pair], pair) for pair in pairs)
    return full_split(best_pair, range(n))

def update_nj(D, index_set):
    """
    Update the distance matrix according to the neighbor joining criterion.
    @param D: the distance matrix
    @param index_set: the subset of indices that will be removed from the updated distance matrix
    @return: an updated distance matrix
    """
    # the neighbor joining update is defined only when two indices are removed
    if len(index_set) != 2:
        raise ValueError('neighbor joining must separate exactly two indices from the remaining indices')
    min_index = min(index_set)
    # define the map from new indices to old indices, with None for the outgroup
    n_big = len(D)
    new_to_old = []
    for i in range(n_big):
        if i == min_index:
            new_to_old.append(None)
        elif i not in index_set:
            new_to_old.append(i)
    # define the outgroup indices
    a_big, b_big = index_set
    # define the new distance matrix
    n_small = len(new_to_old)
    D_small = np.zeros((n_small, n_small))
    for i_small, i_big in enumerate(new_to_old):
        for j_small, j_big in enumerate(new_to_old):
            if i_small == j_small:
                continue
            elif i_big is None:
                D_small[i_small][j_small] = (D[j_big][a_big] + D[j_big][b_big] - D[a_big][b_big])/2
            elif j_big is None:
                D_small[i_small][j_small] = (D[i_big][a_big] + D[i_big][b_big] - D[a_big][b_big])/2
            else:
                D_small[i_small][j_small] = D[i_big][j_big]
    return D_small

def update_generalized_nj(D, index_set):
    """
    Create a new distance matrix according to a neighbor-joining-like criterion.
    Do this according to the explanation in our tree reconstruction manuscript.
    The length of the branch defined by the split is divided evenly
    between the two successor distance matrices.
    @param D: the distance matrix
    @param index_set: the subset of indices that will be removed from the updated distance matrix
    @return: a new distance matrix
    """
    n = len(D)
    A = set(range(n)) - set(index_set)
    B = set(index_set)
    nA = len(A)
    nB = len(B)
    if nA < 2 or nB < 2:
        raise ValueError('expected each side of the split to have at least two elements')
    # The split of the indices into A and B defines a single internal branch.
    # The average distance from A to the branch is alpha.
    # The average distance from B to the branch is beta.
    # The length of the branch is gamma.
    # The expected distance from i to a taxon in the other group is R[i].
    R = {}
    R.update((i, sum(D[i,b] for b in B)/float(nB)) for i in A)
    R.update((j, sum(D[a,j] for a in A)/float(nA)) for j in B)
    gamma_plus_beta = 0.5 * min(R[i]+R[j]-D[i,j] for i, j in itertools.combinations(A, 2))
    alpha_plus_gamma = 0.5 * min(R[i]+R[j]-D[i,j] for i, j in itertools.combinations(B, 2))
    alpha_plus_gamma_plus_beta = sum(D[i,j] for i, j in itertools.product(A, B)) / float(nA * nB)
    gamma = alpha_plus_gamma + gamma_plus_beta - alpha_plus_gamma_plus_beta
    beta = gamma_plus_beta - gamma
    # Initialize the new distance matrix.
    D_out = SchurAlgebra.mmerge(D, index_set)
    # Find the index of D_out that corresponds to the outgroup.
    outgroup_index = sum(1 for a in A if a < min(B))
    D_out[outgroup_index, outgroup_index] = 0
    # Adjust one of the rows and columns to reflect distances to the outgroup.
    label_sets = SchurAlgebra.vmerge([set([i]) for i in range(n)], index_set)
    for i, labels in enumerate(label_sets):
        if i != outgroup_index:
            a = iterutils.get_only(labels)
            d = R[a] - beta - 0.5 * gamma
            D_out[i,outgroup_index] = D_out[outgroup_index,i] = d
    return D_out

def laplacian_to_fiedler(L):
    """
    @param L: the Laplacian matrix
    @return: the Fiedler vector of a related graph
    """
    # get the eigendecomposition
    eigenvalues, V_T = np.linalg.eigh(L)
    eigenvectors = V_T.T.tolist()
    # sort the eigenvectors by their associated eigenvalues
    eigensystem = list(sorted(zip(eigenvalues, eigenvectors)))
    # we are interested in the eigenvector whose eigenvalue is second least
    w, v = eigensystem[1]
    return v

def dccov_to_fiedler(HSH):
    """
    @param HSH: the doubly centered covariance matrix
    @return: the Fiedler vector of a related graph
    """
    # get the eigendecomposition
    eigenvalues, V_T = np.linalg.eigh(HSH)
    eigenvectors = V_T.T.tolist()
    # get the eigenvalue and eigenvector of interest
    w, v = max(zip(eigenvalues, eigenvectors))
    return v

def edm_to_fiedler(D):
    """
    @param D: the distance matrix
    @return: the Fiedler vector of a related graph
    """
    return dccov_to_fiedler(Euclid.edm_to_dccov(D))

def split_using_eigenvector(D, epsilon=1e-14):
    """
    Split the distance matrix using signs of an eigenvector of -HDH/2.
    If a degenerate split is found then a DegenerateSplitException is raised.
    @param D: the distance matrix
    @param epsilon: small eigenvector loadings will be treated as zero
    @return: a set of two index sets defining a split of the indices
    """
    # get the fiedler vector
    v = edm_to_fiedler(D)
    # get the eigensplit
    eigensplit = eigenvector_to_split(v, epsilon)
    # validate the split
    min_cardinality, min_set = min((len(s), s) for s in eigensplit)
    if min_cardinality == 0:
        raise InvalidSpectralSplitException(D)
    elif min_cardinality == 1:
        index, = min_set
        raise DegenerateSplitException(index)
    else:
        return eigensplit

def split_using_eigenvector_with_nj_fallback(D):
    """
    Try to split using an eigenvector, but fall back to a neighbor joining split if the split is degenerate.
    @param D: the distance matrix
    @return: a set of two index sets defining a split of the indices
    """
    try:
        return split_using_eigenvector(D)
    except DegenerateSplitException, e:
        return split_nj(D)

def update_using_laplacian(D, index_set):
    """
    Update the distance matrix by summing rows and columns of the removed indices.
    @param D: the distance matrix
    @param index_set: the set of indices that will be removed from the updated distance matrix
    @return: an updated distance matrix
    """
    L = Euclid.edm_to_laplacian(D)
    L_small = SchurAlgebra.mmerge(L, index_set)
    D_small = Euclid.laplacian_to_edm(L_small)
    return D_small


class TestBuildTreeTopology(unittest.TestCase):

    def setUp(self):
        """
        Define a perturbed and a true distance matrix from a paper.
        The paper is Why Neighbor Joining Works.
        """
        self.D_perturbed = np.array([
            [  0, 2.7, 2.6, 2.6, 2.6, 4.4, 4.4, 4.4],
            [2.7,   0, 4.4, 4.4, 4.4, 2.6, 2.6, 2.6],
            [2.6, 4.4,   0, 0.1, 0.4, 2.7, 2.7, 2.7],
            [2.6, 4.4, 0.1,   0, 0.4, 2.7, 2.7, 2.7],
            [2.6, 4.4, 0.4, 0.4,   0, 2.7, 2.7, 2.7],
            [4.4, 2.6, 2.7, 2.7, 2.7,   0, 0.1, 0.4],
            [4.4, 2.6, 2.7, 2.7, 2.7, 0.1,   0, 0.4],
            [4.4, 2.6, 2.7, 2.7, 2.7, 0.4, 0.4,   0]])
        self.D = np.array([
            [0.0, 3.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0],
            [3.0, 0.0, 3.0, 3.0, 3.0, 2.0, 2.0, 2.0],
            [2.0, 3.0, 0.0, 0.1, 0.4, 3.0, 3.0, 3.0],
            [2.0, 3.0, 0.1, 0.0, 0.4, 3.0, 3.0, 3.0],
            [2.0, 3.0, 0.4, 0.4, 0.0, 3.0, 3.0, 3.0],
            [3.0, 2.0, 3.0, 3.0, 3.0, 0.0, 0.1, 0.4],
            [3.0, 2.0, 3.0, 3.0, 3.0, 0.1, 0.0, 0.4],
            [3.0, 2.0, 3.0, 3.0, 3.0, 0.4, 0.4, 0.0]])
        abc = (4, 0.2), (((2, 0.05), (3, 0.05)), 0.15)
        mnp = (7, 0.2), (((5, 0.05), (6, 0.05)), 0.15)
        self.tree = Xtree.list_to_weighted_tree(((((abc, 0.8), (0, 1.0)), 1.0), (1, 1.0), (mnp, 0.8)))
        self.true_splits = set([
            make_split((2,3),(0,1,4,5,6,7)),
            make_split((2,3,4),(0,1,5,6,7)),
            make_split((2,3,4,0),(1,5,6,7)),
            make_split((2,3,4,0,1),(5,6,7)),
            make_split((2,3,4,0,1,7),(5,6))])

    def test_true_matrix(self):
        """
        Each method should work on the true distance matrix.
        """
        # assert that neighbor joining finds the correct splits
        observed = get_splits(self.D, split_nj, update_nj)
        self.assertEqual(observed, self.true_splits)
        # assert that spectral splitting with neighbor joining fallback finds the correct splits
        observed = get_splits(self.D, split_using_eigenvector_with_nj_fallback, update_using_laplacian)
        self.assertEqual(observed, self.true_splits)

    def test_perturbed_matrix(self):
        """
        Only the spectral method should work on the perturbed distance matrix.
        """
        # assert that neighbor joining fails to find the correct splits of the perturbed matrix
        observed = get_splits(self.D_perturbed, split_nj, update_nj)
        self.assertNotEqual(observed, self.true_splits)
        # assert that spectral splitting with neighbor joining fallback finds the correct splits
        observed = get_splits(self.D_perturbed, split_using_eigenvector_with_nj_fallback, update_using_laplacian)
        self.assertEqual(observed, self.true_splits)

    def test_balanced_tree(self):
        """
        Test splits of a balanced tree.
        The tree from Why Neighbor Joining Works is a caterpillar tree.
        """
        D_balanced = np.array([
            [0, 2, 4, 4, 5, 5, 5, 5],
            [2, 0, 4, 4, 5, 5, 5, 5],
            [4, 4, 0, 2, 5, 5, 5, 5],
            [4, 4, 2, 0, 5, 5, 5, 5],
            [5, 5, 5, 5, 0, 2, 4, 4],
            [5, 5, 5, 5, 2, 0, 4, 4],
            [5, 5, 5, 5, 4, 4, 0, 2],
            [5, 5, 5, 5, 4, 4, 2, 0]])
        true_splits = set([
            make_split((0,1),(2,3,4,5,6,7)),
            make_split((2,3),(0,1,4,5,6,7)),
            make_split((4,5),(0,1,2,3,6,7)),
            make_split((6,7),(0,1,2,3,4,5)),
            make_split((0,1,2,3),(4,5,6,7))])
        # assert that neighbor joining finds the correct splits
        observed = get_splits(D_balanced, split_nj, update_nj)
        self.assertEqual(observed, true_splits)
        # assert that spectral splitting with neighbor joining fallback finds the correct splits
        observed = get_splits(D_balanced, split_using_eigenvector_with_nj_fallback, update_using_laplacian)
        self.assertEqual(observed, true_splits)

    def test_quartet_additivity(self):
        """
        Test quartet additivity.
        """
        self.assertFalse(is_quartet_additive(self.tree, self.D_perturbed))
        self.assertTrue(is_quartet_additive(self.tree, self.D))

    def test_quartet_consistency(self):
        """
        Test quartet additivity.
        """
        self.assertTrue(is_quartet_consistent(self.tree, self.D))

    def test_update_generalized_nj_small(self):
        """
        Test the generation of successor distance matrices from a small initial distance matrix.
        """
        D = np.array([
            [0, 3, 8, 9],
            [3, 0, 9, 10],
            [8, 9, 0, 9],
            [9, 10, 9, 0]])
        # check the first successor distance matrix
        D_out_a_expected = ([
            [0, 3, 2.5],
            [3, 0, 3.5],
            [2.5, 3.5, 0]])
        index_set = set([2,3])
        D_out_a = update_generalized_nj(D, index_set)
        self.assertTrue(np.allclose(D_out_a_expected, D_out_a), msg=D_out_a)
        # check the second successor distance matrix
        D_out_b_expected = ([
            [0, 5.5, 6.5],
            [5.5, 0, 9],
            [6.5, 9, 0]])
        index_set = set([0,1])
        D_out_b = update_generalized_nj(D, index_set)
        self.assertTrue(np.allclose(D_out_b_expected, D_out_b), msg=D_out_b)

    def test_update_generalized_nj_big(self):
        """
        Test the generation of successor distance matrices from a more complicated initial distance matrix.
        """
        # define the initial tree and the two subtrees
        s_tree_initial = '(((3:9, 2:2):4, 1:2):1, (4:1, 5:3):7, 6:2);'
        s_tree_a = '((3:9, 2:2):4, 1:2, B:0.5);'
        s_tree_b = '((4:1, 5:3):7, 6:2, A:0.5);'
        # Define an ordering of the taxa.
        # The initial ordering is arbitrary,
        # and the subsequent orderings are dependent on the initial ordering.
        taxa_initial = ['1', '4', '2', '5', '3', '6']
        taxa_a = ['1', 'B', '2', '3']
        taxa_b = ['A', '4', '5', '6']
        # Define the distance matrices.
        D_initial = np.array(NewickIO.parse(s_tree_initial, FelTree.NewickTree).get_distance_matrix(taxa_initial))
        D_a = np.array(NewickIO.parse(s_tree_a, FelTree.NewickTree).get_distance_matrix(taxa_a))
        D_b = np.array(NewickIO.parse(s_tree_b, FelTree.NewickTree).get_distance_matrix(taxa_b))
        # assert that the correct distance matrices are created
        D_out_a = update_generalized_nj(D_initial, set([1,3,5]))
        D_out_b = update_generalized_nj(D_initial, set([0,2,4]))
        self.assertTrue(np.allclose(D_a, D_out_a))
        self.assertTrue(np.allclose(D_b, D_out_b))


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestBuildTreeTopology)
    unittest.TextTestRunner(verbosity=2).run(suite)
