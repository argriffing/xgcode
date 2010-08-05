"""
Implmentation of the gap statistic for clustering.

This is implemented from the description
in the paper entitled 'estimating the number of clusters
in a data set via the gap statistic' by Tibshirani et al.
"""

import unittest
import itertools
import random

import numpy as np

def points_to_ranges(X):
    """
    @param X: each row of numpy array X is a point
    @return: a list giving, for each coordinate, a pair of low, high bounds
    """
    lows = np.min(X, axis=0)
    highs = np.max(X, axis=0)
    return zip(lows.tolist(), highs.tolist())

def sample_points(ranges, n):
    """
    @param ranges: for each coordinate, a pair of low, high bounds
    @param n: the number of points to sample
    @return: a 2D numpy array where each row is a point
    """
    return np.array([_sample_point(ranges) for i in range(n)])

def get_simulation_correction(wlogs):
    """
    The input is the array of logs of wcss values from the simulations.
    This is defined in the third step of the computation of the
    gap statistic in section four of the gap statistic paper.
    @param wlogs: numpy array of logs of within cluster sums of squares
    @return: the correction for simulation error
    """
    B = len(ws)
    sd = np.std(wlogs)
    return sd * math.sqrt(1.0 + 1.0/B)

def get_wcss(X):
    """
    @param X: each row of numpy array X is a point
    @return: sum of squares of distances to the center
    """
    nrows, ncols = X.shape
    center = X.mean(axis=0)
    return np.sum((X-center)**2)

def get_sum_of_distances(X):
    """
    @param X: each row of numpy array X is a point
    @return: sum of squares of euclidean distances between all ordered pairs
    """
    return _get_sum_of_distances_fast(X)

def _sample_point(ranges):
    """
    @param ranges: for each coordinate, a pair of low, high bounds
    @return: a float sequence
    """
    return [random.uniform(low, high) for low, high in ranges]

def _get_sum_of_distances_fast(X):
    """
    @param X: each row of numpy array X is a point
    @return: sum of squares of euclidean distances between all ordered pairs
    """
    nrows, ncols = X.shape
    colsums = X.sum(axis=0)
    alpha = nrows * np.sum(X**2)
    beta = np.dot(colsums, colsums)
    return 2*(alpha - beta)

def _get_sum_of_distances_centered(X):
    """
    @param X: each row of numpy array X is a point
    @return: sum of squares of euclidean distances between all ordered pairs
    """
    nrows, ncols = X.shape
    center = X.mean(axis=0)
    return 2*nrows*np.sum((X-center)**2)

def _get_sum_of_distances_naive(X):
    """
    @param X: each row of numpy array X is a point
    @return: sum of squares of euclidean distances between all ordered pairs
    """
    return sum(np.dot(a-b, a-b) for a, b in itertools.product(X, X))


class TestTibshirani(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        unittest.TestCase.__init__(self, *args, **kwargs)
        self.M = np.array([
            [1, 2, 3, 4],
            [1, 1, 1, 1],
            [8, 2, 3, 6]], dtype=float)

    def assertAllclose(self, x, y):
        self.assert_(np.allclose(x, y), (x, y)) 

    def test_sum_of_distances(self):
        points = np.array([[0,0], [3,4]], dtype=float)
        observed = get_sum_of_distances(points)
        self.assertAllclose(observed, 50.0)

    def test_sum_of_distances_fast(self):
        expected = _get_sum_of_distances_naive(self.M)
        observed = _get_sum_of_distances_fast(self.M)
        self.assertTrue(np.allclose(expected, observed))
    
    def test_sum_of_distances_centered(self):
        expected = _get_sum_of_distances_naive(self.M)
        observed = _get_sum_of_distances_centered(self.M)
        self.assertTrue(np.allclose(expected, observed))

if __name__ == '__main__':
    unittest.main()
