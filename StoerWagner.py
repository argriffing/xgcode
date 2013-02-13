"""
A polynomial min cut algorithm.
"""

from StringIO import StringIO
import unittest

import numpy as np

import MatrixUtil


# This matrix is an example used in the Stoer-Wagner min cut paper.
g_stoer_wagner_affinity = np.array([
        [0, 2, 0, 0, 3, 0, 0, 0],
        [2, 0, 3, 0, 2, 2, 0, 0],
        [0, 3, 0, 4, 0, 0, 2, 0],
        [0, 0, 4, 0, 0, 0, 2, 2],
        [3, 2, 0, 0, 0, 3, 0, 0],
        [0, 2, 0, 0, 3, 0, 1, 0],
        [0, 0, 2, 2, 0, 1, 0, 3],
        [0, 0, 0, 2, 0, 0, 3, 0],
        ], dtype=float)

def stoer_wagner_min_cut(weight_matrix):
    """
    The input matrix is assumed to be a numpy ndarray with float dtype.
    @param weight_matrix: non-negative symmetric weighted adjacency matrix
    @return: the set of indices belonging to one of the two clusters
    """
    w = weight_matrix.copy()
    n = w.shape[0]
    MatrixUtil.assert_symmetric(w)
    MatrixUtil.assert_nonnegative(w)
    MatrixUtil.assert_hollow(w)
    # no cut has been observed so far
    min_cut = None
    min_cut_weight = None
    remaining_indices = set(range(n))
    index_sets = [set([i]) for i in range(n)]
    # reduce the number of remaining indices by one each iteration
    while len(remaining_indices) > 1:
        # initialize the borg
        assimilated_indices = [0]
        unassimilated_indices = remaining_indices - set(assimilated_indices)
        weight_index_pairs = [(w[i, 0], i) for i in unassimilated_indices]
        while len(assimilated_indices) < len(remaining_indices):
            max_weight, best_index = max(weight_index_pairs)
            assimilated_indices.append(best_index)
            unassimilated_indices.remove(best_index)
            next_weight_index_pairs = []
            for weight, index in weight_index_pairs:
                if index is not best_index:
                    next_weight = weight + w[index, best_index]
                    next_weight_index_pairs.append((next_weight, index))
            weight_index_pairs = next_weight_index_pairs
        # the cut between the last two assimilated indices is a possible min cut
        cut = set(index_sets[assimilated_indices[-1]])
        s, t = list(sorted(assimilated_indices[-2:]))
        cut_weight = max_weight
        if min_cut is None or cut_weight < min_cut_weight:
            min_cut, min_cut_weight = cut, cut_weight
        # combine the last two assimilated indices
        for i in remaining_indices:
            w[i, s] += w[i, t]
            w[s, i] += w[t, i]
        index_sets[s].update(index_sets[t])
        remaining_indices.remove(t)
    return min_cut


class TestStoerWagner(unittest.TestCase):

    def test_stoer_wagner_example(self):
        w = g_stoer_wagner_affinity
        n = len(w)
        cut = stoer_wagner_min_cut(w)
        expected_cuts = (set([2, 3, 6, 7]), set([0, 1, 4, 5]))
        if cut not in expected_cuts:
            self.fail('the cut was not correct')
        weight = 0
        for i in cut:
            for j in set(range(n)) - cut:
                weight += w[i, j]
        self.assertEqual(weight, 4)


if __name__ == '__main__':
    unittest.main()
