"""
A polynomial min cut algorithm.
"""

from StringIO import StringIO
import unittest


# This matrix is an example used in the Stoer-Wagner min cut paper.
g_stoer_wagner_affinity = [
        [0, 2, 0, 0, 3, 0, 0, 0],
        [2, 0, 3, 0, 2, 2, 0, 0],
        [0, 3, 0, 4, 0, 0, 2, 0],
        [0, 0, 4, 0, 0, 0, 2, 2],
        [3, 2, 0, 0, 0, 3, 0, 0],
        [0, 2, 0, 0, 3, 0, 1, 0],
        [0, 0, 2, 2, 0, 1, 0, 3],
        [0, 0, 0, 2, 0, 0, 3, 0]]

def stoer_wagner_min_cut(weight_matrix):
    """
    @param weight_matrix: a lists of lists that defines a symmetric matrix
    @return: the set of indices belonging to one of the two clusters
    """
    # create a copy of the weight matrix so it can be modified
    w = [row[:] for row in weight_matrix]
    n = len(w)
    # validate the weight matrix
    for i in range(n):
        for j in range(n):
            if weight_matrix[i][j] < 0:
                raise ValueError('each element of the weight matrix must be non-negative')
    for i in range(n):
        if weight_matrix[i][i] != 0:
            raise ValueError('each diagonal element of the weight matrix must be zero')
    for i in range(n):
        for j in range(i):
            if weight_matrix[i][j] != weight_matrix[j][i]:
                raise ValueError('the weight matrix must be symmetric')
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
        weight_index_pairs = [(w[i][0], i) for i in unassimilated_indices]
        while len(assimilated_indices) < len(remaining_indices):
            max_weight, best_index = max(weight_index_pairs)
            assimilated_indices.append(best_index)
            unassimilated_indices.remove(best_index)
            next_weight_index_pairs = []
            for weight, index in weight_index_pairs:
                if index is not best_index:
                    next_weight = weight + w[index][best_index]
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
            w[i][s] += w[i][t]
            w[s][i] += w[t][i]
        index_sets[s].update(index_sets[t])
        remaining_indices.remove(t)
    return min_cut


class TestStoerWagner(unittest.TestCase):

    def test_stoer_wagner_example(self):
        matrix = g_stoer_wagner_affinity
        n = len(matrix)
        cut = stoer_wagner_min_cut(matrix)
        expected_cuts = (set([2, 3, 6, 7]), set([0, 1, 4, 5]))
        if cut not in expected_cuts:
            self.fail('the cut was not correct')
        weight = 0
        for i in cut:
            for j in set(range(n)) - cut:
                weight += matrix[i][j]
        self.assertEqual(weight, 4)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestStoerWagner)
    unittest.TextTestRunner(verbosity=2).run(suite)
