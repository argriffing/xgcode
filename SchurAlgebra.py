"""
This module defines analogous functions on square matrices and vectors of disjoint integer sets.
"""

import unittest

import numpy as np

def assert_ordered(label_sets):
    """
    Assert an invariant of the list of label sets.
    @param label_sets: an ordered list of sets of disjoint integers
    """
    min_elements = [min(labels) for labels in label_sets]
    for a, b in zip(min_elements[:-1], min_elements[1:]):
        assert a < b

def assert_disjoint(label_sets):
    """
    Assert an invariant of the list of label sets.
    @param label_sets: an ordered list of sets of disjoint integers
    """
    # get the total number of elements
    n = sum(len(s) for s in label_sets)
    # get the number of unique elements
    total_union = set()
    for s in label_sets:
        total_union.update(s)
    # assert that each element is unique
    assert n == len(total_union)

def assert_valid_indices(n, index_set):
    """
    Assert that all indices are in bounds.
    @param n: valid indices are in range(n)
    @param index_set: a set of valid indices
    """
    assert index_set <= set(range(n))

def assert_square(M):
    """
    Assert that a matrix is square.
    @param M: a matrix
    """
    s = M.shape
    if len(s) != 2 or s[0] != s[1]:                                
        assert False, 'the matrix is not square'

def schur_helper(M, nremove):
    """
    Remove some number of final rows and columns by schur complementation.
    @param M: a square rank 2 numpy array
    @param nremove: Schur out this many of the last rows and columns
    @return: a new square rank 2 numpy array
    """
    assert_square(M)
    if not nremove:
        return M
    n = len(M)
    nkeep = n - nremove
    # If the matrix is symmetric then B should be transpose of D.
    A = M[:-nremove][:, :-nremove]
    B = M[:-nremove][:, nkeep:]
    C = M[nkeep:][:, nkeep:]
    D = M[nkeep:][:, :-nremove]
    C_pinv = np.linalg.pinv(C)
    return A - np.dot(B, np.dot(C_pinv, D))

def vdelete(label_sets, index_set):
    """
    @param label_sets: an ordered list of sets of disjoint integers
    @param index_set: a set of indices to delete
    @return: a stably reduced list of label sets
    """
    assert_ordered(label_sets)
    assert_disjoint(label_sets)
    assert_valid_indices(len(label_sets), index_set)
    return [s for i, s in enumerate(label_sets) if i not in index_set]

def vmerge(label_sets, index_set):
    """
    @param label_sets: an ordered list of sets of disjoint integers
    @param index_set: a set of indices to merge
    @return: a stably merged list of label sets
    """
    assert_ordered(label_sets)
    assert_disjoint(label_sets)
    assert_valid_indices(len(label_sets), index_set)
    # define the merged set of labels
    merged_labels = set()
    for i in index_set:
        merged_labels.update(label_sets[i])
    # get the min index which will be replaced
    min_index = min(index_set)
    # define the new list of label sets
    next_label_sets = []
    for i, label_set in enumerate(label_sets):
        if i == min_index:
            next_label_sets.append(merged_labels)
        elif i not in index_set:
            next_label_sets.append(label_set)
    return next_label_sets

def mdelete(M, index_set):
    """
    @param M: a symmetric square matrix
    @param index_set: a set of indices to delete
    @return: a symmetric square matrix with some rows and columns deleted
    """
    assert_square(M)
    assert_valid_indices(len(M), index_set)
    ind = [i for i in range(len(M)) if i not in index_set]
    return M[np.ix_(ind, ind)].copy()

def mmerge(M, index_set):
    """
    @param M: a symmetric square matrix
    @param index_set: a set of indices to merge
    """
    assert_square(M)
    assert_valid_indices(len(M), index_set)
    n = len(M)
    index_sets = vmerge([set([i]) for i in range(n)], index_set)
    n_small = len(index_sets)
    M_small = np.zeros((n_small, n_small))
    for i, set_i in enumerate(index_sets):
        for j, set_j in enumerate(index_sets):
            M_small[i][j] = sum(M[k][l] for k in set_i for l in set_j)
    return M_small

def mschur(M, index_set):
    """
    @param M: a symmetric square matrix
    @param index_set: a set of indices to be removed by schur complementation
    """
    assert_square(M)
    assert_valid_indices(len(M), index_set)
    n = len(M)
    total_set = set(range(n))
    remove_set = set(index_set)
    keep_set = total_set - remove_set
    # define the permutation
    new_to_old = list(sorted(keep_set)) + list(sorted(remove_set))
    # get the permuted matrix
    M_permuted = np.zeros_like(M)
    for new_i, old_i in enumerate(new_to_old):
        for new_j, old_j in enumerate(new_to_old):
            M_permuted[new_i][new_j] = M[old_i][old_j]
    # schur complement out the last rows and columns
    return schur_helper(M_permuted, len(remove_set))


class TestSchurAlgebra(unittest.TestCase):

    def test_vdelete(self):
        label_sets = [set([1, 2, 3]), set([4, 5]), set([6]), set([7, 8])]
        index_set = set([1, 2])
        expected = [set([1, 2, 3]), set([7, 8])]
        observed = vdelete(label_sets, index_set)
        self.assertEquals(observed, expected)

    def test_vmerge(self):
        # setup
        label_sets = [set([1, 2, 3]), set([4, 5]), set([6]), set([7, 8])]
        # first test
        index_set = set([1, 2])
        expected = [set([1, 2, 3]), set([4, 5, 6]), set([7, 8])]
        observed = vmerge(label_sets, index_set)
        self.assertEquals(observed, expected)
        # second test
        index_set = set([2, 3])
        expected = [set([1, 2, 3]), set([4, 5]), set([6, 7, 8])]
        observed = vmerge(label_sets, index_set)
        self.assertEquals(observed, expected)
        # third test
        index_set = set([0, 1])
        expected = [set([1, 2, 3, 4, 5]), set([6]), set([7, 8])]
        observed = vmerge(label_sets, index_set)
        self.assertEquals(observed, expected)

    def test_mmerge(self):
        M = np.array([
            [0, 1, 2, 3],
            [4, 5, 6, 7],
            [8, 9, 0, 1],
            [2, 3, 4, 5]])
        # first test
        index_set = set([0, 1])
        expected = np.array([
            [0+1+4+5, 2+6, 3+7],
            [8+9, 0, 1],
            [2+3, 4, 5]])
        observed = mmerge(M, index_set)
        self.assertTrue(np.allclose(observed, expected))

    def test_mschur(self):
        """
        Some of these tests are taken from a book by Fuzhen Zhang.
        The title is "The Schur complement and its applications".
        """
        # a test from section 4.2
        M = np.array([
            [1, 1, 1],
            [-1, 1, 2],
            [0, -1, 1]])
        index_set = set([1])
        expected = np.array([
            [2, -1],
            [-1, 3]])
        observed = mschur(M, index_set)
        self.assertTrue(np.allclose(observed, expected))
        # a test from section 4.2
        M = np.array([
            [0, 1, 2, 2],
            [1, 0, 1, 5],
            [2, 1, 0, 10],
            [2, 5, 10, 0]])
        index_set = set([0, 1])
        expected = np.array([
            [-4, -2],
            [-2, -20]])
        observed = mschur(M, index_set)
        self.assertTrue(np.allclose(observed, expected))

    def test_mdelete(self):
        M = np.array([
            [0, 1, 2, 3],
            [4, 5, 6, 7],
            [8, 9, 0, 1],
            [2, 3, 4, 5]])
        index_set = set([1, 2])
        expected = np.array([
            [0, 3],
            [2, 5]])
        observed = mdelete(M, index_set)
        self.assertTrue(np.allclose(observed, expected))

if __name__ == '__main__':
    unittest.main()
