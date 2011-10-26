"""
This module uses divisive clustering without outgrouping to define a tree topology.
This is a naive method to be contrasted with the more advanced outgrouping approach.
"""

import unittest
import random

import numpy as np

import SchurAlgebra
import Euclid
import BuildTreeTopology
import MatrixUtil

def random_split(D):
    """
    @param D: numpy distance matrix
    @return: a split of the indices of the distance matrix
    """
    n = len(D)
    nsample = random.randrange(1, n)
    left = frozenset(random.sample(range(n), nsample))
    right = frozenset(range(n)) - left
    return frozenset([left, right])

def spectral_split(D):
    """
    @param D: numpy distance matrix
    @return: a split of the indices of the distance matrix
    """
    v = BuildTreeTopology.edm_to_fiedler(D)
    split = BuildTreeTopology.eigenvector_to_split(v)
    # if one side of the split is empty then there is a failure
    if frozenset() in split:
        raise RuntimeError('failed split')
    return split

def get_hierarchy(D, split_function, labels, nlevel_limit=None):
    """
    @param D: numpy distance matrix
    @param split_function: returns a label split given a distance matrix and a list of labels
    @param labels: an ordered list of labels
    @param nlevel_limit: the depth to which the indices should be split
    @return: a nested frozenset of indices of the original distance matrix
    """
    # sanity check
    assert len(D) == len(labels)
    # if labels cannot be split further then stop
    if len(labels) == 1:
        return labels[0]
    # if we have reached the requested depth then stop
    if nlevel_limit == 0:
        return frozenset(labels)
    # if there are only two labels then there can be only one split
    if len(labels) == 2:
        alpha, beta = labels
        return frozenset([frozenset([alpha]), frozenset([beta])])
    # get the next split limit
    next_nlevel_limit = None if nlevel_limit is None else nlevel_limit - 1
    # do the split
    index_split = split_function(D)
    # get each side of the split recursively
    sides = []
    for index_set in index_split:
        indices = list(sorted(index_set))
        D_next = MatrixUtil.get_principal_submatrix(D, indices)
        labels_next = [labels[i] for i in indices]
        sides.append(get_hierarchy(D_next, split_function, labels_next, next_nlevel_limit))
    return frozenset(sides)

def _build_clusters(hierarchy, clusters):
    """
    This is a helper function.
    @param hierarchy: a nested frozenset of indices
    @param clusters: the current list of clusters
    @return: the flattened hierarchy
    """
    all_indices = set()
    try:
        for sub_element in hierarchy:
            all_indices.update(_build_clusters(sub_element, clusters))
    except TypeError, e:
        return set([hierarchy])
    clusters.append(all_indices)
    return all_indices

def hierarchy_to_nontrivial_splits(hierarchy):
    clusters = []
    all_indices = _build_clusters(hierarchy, clusters)
    for cluster in clusters:
        left = frozenset(cluster)
        right = frozenset(all_indices - cluster)
        split = frozenset([left, right])
        if min(len(side) for side in split) > 1:
            yield split


class TestDendrogram(unittest.TestCase):

    def test_3x3(self):
        """
        This tree is constructed correctly.
        """
        D = np.array([
            [0, 3, 3],
            [3, 0, 2],
            [3, 2, 0]])
        labels = [100, 200, 300]
        observed_hierarchy = get_hierarchy(D, spectral_split, labels)
        observed_nontrivial_splits = set(hierarchy_to_nontrivial_splits(observed_hierarchy))
        # check the hierarcy
        expected_hierarchy = frozenset([100, frozenset([200, 300])])
        self.assertEqual(expected_hierarchy, observed_hierarchy)
        # check the nontrivial splits
        expected_nontrivial_splits = set()
        self.assertEqual(expected_nontrivial_splits, observed_nontrivial_splits)

    def test_4x4(self):
        """
        This tree is constructed correctly.
        """
        D = np.array([
            [0, 2, 3, 3],
            [2, 0, 3, 3],
            [3, 3, 0, 2],
            [3, 3, 2, 0]])
        labels = [100, 200, 300, 400]
        observed_hierarchy = get_hierarchy(D, spectral_split, labels)
        observed_nontrivial_splits = set(hierarchy_to_nontrivial_splits(observed_hierarchy))
        # check the hierarcy
        expected_hierarchy = frozenset([frozenset([100, 200]), frozenset([300, 400])])
        self.assertEqual(expected_hierarchy, observed_hierarchy)
        # check the nontrivial splits
        expected_nontrivial_splits = set([frozenset([frozenset([100, 200]), frozenset([300, 400])])])
        self.assertEqual(expected_nontrivial_splits, observed_nontrivial_splits)

    def test_5x5(self):
        """
        This tree is constructed incorrectly.
        """
        # true tree:
        # ((a:1, b:2):1, c:1, (d:1, e:1):10);
        # naively reconstructed topology:
        # ((a, c), b, (d, e))
        D = np.array([
            [0, 6, 3, 13, 13],
            [6, 0, 7, 17, 17],
            [3, 7, 0, 12, 12],
            [13, 17, 12, 0, 2],
            [13, 17, 12, 2, 0]])
        labels = [100, 200, 300, 400, 500]
        observed_hierarchy = get_hierarchy(D, spectral_split, labels)
        observed_nontrivial_splits = set(hierarchy_to_nontrivial_splits(observed_hierarchy))
        # check the hierarcy
        expected_hierarchy = frozenset([frozenset([200, frozenset([300, 100])]), frozenset([400, 500])])
        self.assertEqual(expected_hierarchy, observed_hierarchy)
        # check the nontrivial splits
        expected_nontrivial_splits = set([
            frozenset([frozenset([200, 100, 300]), frozenset([400, 500])]),
            frozenset([frozenset([200, 400, 500]), frozenset([300, 100])])])
        self.assertEqual(expected_nontrivial_splits, observed_nontrivial_splits)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestDendrogram)
    unittest.TextTestRunner(verbosity=2).run(suite)
