"""
Hierarchical agglomerative clustering.
The main data structures are a list of points, a map from a
cluster index to a set of point indices, and a map from an ordered (low,high)
pair of cluster indices to a sum of squared distances.
The motivation for this module is to do clustering in a simple deterministic
way so that an extra layer of cluster significance testing can be added later.
Therefore it uses squared error and average linkage, following
example 4.1 in the gap statistic paper by Tibshirani et al.
The terms ssd and msd are used for sum of squared Euclidean distances
and mean of squared Euclidean distances respectively.
"""

import unittest
import itertools
import sys

import numpy as np


def get_initial_cluster_map(points):
    """
    @param points: a sequence of points in euclidean space
    @return: an initial cluster map
    """
    return dict((i, set([i])) for i, p in enumerate(points))

def get_initial_b_ssd_map(points):
    """
    @param points: a sequence of points in euclidean space
    @return: an initial between-clusters ssd map
    """
    b_ssd_map = {}
    npoints = len(points)
    for j in range(npoints):
        for i in range(j):
            diff = points[j] - points[i]
            b_ssd_map[i, j] = np.dot(diff, diff)
    return b_ssd_map

def get_initial_w_ssd_map(points):
    """
    @param points: a sequence of points in euclidean space
    @return: an initial within-clusters ssd map
    """
    return dict((i, 0.0) for i, p in enumerate(points))

def merge(cluster_map, w_ssd_map, b_ssd_map, pair):
    """
    Modify the maps in-place.
    The new cluster index is one greater than the max existing cluster index.
    @param cluster_map: maps a cluster index to a set of point indices
    @param w_ssd_map: maps a cluster index to a sum of squared distances
    @param b_ssd_map: maps a cluster index pair to a sum of squared distances
    @param pair: the pair of cluster indices to merge
    """
    # define the cluster indices to be merged
    ia, ib = pair
    # define the index of the merged cluster
    neo = max(cluster_map) + 1
    # define the ssd keys that contain a constituent cluster index
    dead_ssd_keys = set((i, j) for i, j in b_ssd_map if set([i,j]) & set(pair))
    # add ssd to the new cluster to the between cluster sum of squares
    for i in cluster_map:
        if i in pair:
            # do not compute distances to the constituent clusters
            continue
        ka = (i,ia) if i < ia else (ia,i)
        kb = (i,ib) if i < ib else (ib,i)
        b_ssd_map[i,neo] = b_ssd_map[ka] + b_ssd_map[kb]
    # add the new cluster to the cluster map
    cluster_map[neo] = cluster_map[ia] | cluster_map[ib]
    # add the new cluster to the within cluster sum of squares
    kw = (ia,ib) if ia < ib else (ib,ia)
    w_ssd_map[neo] = w_ssd_map[ia] + w_ssd_map[ib] + b_ssd_map[kw]
    # remove the constituent clusters from the between ssd map
    for k in dead_ssd_keys:
        del b_ssd_map[k]
    # remove the constituent clusters from the cluster map
    for k in pair:
        del cluster_map[k]
    # remove the constituent clusters from the within ssd map
    for k in pair:
        del w_ssd_map[k]

def get_msd(cluster_map, b_ssd_map, pair):
    """
    This function should take constant time.
    @param cluster_map: maps a cluster index to a set of point indices
    @param b_ssd_map: maps a cluster index pair to a sum of squared distances
    @param pair: the pair of cluster indices for which to compute msd
    """
    i, j = pair
    return b_ssd_map[pair] / float(len(cluster_map[i]) + len(cluster_map[j]))

def get_pair(cluster_map, b_ssd_map):
    """
    Get the indices of the best pair of clusters to merge.
    The best pair is the pair with the lowest mean squared distance.
    @param cluster_map: maps a cluster index to a set of point indices
    @param ssd_map: maps a cluster index pair to a sum of squared distances
    @return: the pair of indices
    """
    msd, p = min((get_msd(cluster_map, b_ssd_map, p), p) for p in b_ssd_map)
    return p


class TestAgglom(unittest.TestCase):

    def test_placeholder(self):
        pass


if __name__ == '__main__':
    unittest.main()
