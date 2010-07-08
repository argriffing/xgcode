"""
Implmentation of the gap statistic for clustering.

This is implemented from the description in the paper entitled
estimating the number of clusters in a data set via the gap statistic.
"""

import unittest
import itertools

import numpy as np

def get_sum_of_distances_fast(X):
    """
    @param X: each row of numpy array X is a point
    @return: sum of euclidean distances between all ordered pairs
    """
    nrows, ncols = X.shape
    colsums = X.sum(axis=0)
    alpha = nrows * np.sum(X**2)
    beta = np.dot(colsums, colsums)
    return 2 * (alpha - beta)

def get_sum_of_distances_naive(X):
    """
    @param X: each row of numpy array X is a point
    @return: sum of euclidean distances between all ordered pairs
    """
    return sum(np.dot(a-b, a-b) for a, b in itertools.product(X, X))


class TestTibshirani(unittest.TestCase):
    
    def test_sum_of_distances(self):
        M = np.array([
            [1, 2, 3, 4],
            [1, 1, 1, 1],
            [8, 2, 3, 6]])
        expected = get_sum_of_distances_naive(M)
        observed = get_sum_of_distances_fast(M)
        self.assertTrue(np.array_equal(expected, observed))

