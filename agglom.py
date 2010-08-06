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
import heapq

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
    dead_ssd_keys = get_dead_ssd_keys(cluster_map, pair)
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
    @param b_ssd_map: maps a cluster index pair to a sum of squared distances
    @return: the pair of indices
    """
    msd, p = min((get_msd(cluster_map, b_ssd_map, p), p) for p in b_ssd_map)
    return p

def get_dead_ssd_keys(cluster_map, pair):
    """
    @param cluster_map: maps a cluster index to a set of point indices
    @param pair: the pair of merged indices
    """
    ia, ib = pair
    dead_ssd_keys = set()
    for i in cluster_map:
        if i != ia:
            dead_ssd_keys.add(tuple(sorted((i,ia))))
        if i != ib:
            dead_ssd_keys.add(tuple(sorted((i,ib))))
    dead_ssd_keys.add(tuple(sorted(pair)))
    return dead_ssd_keys

def merge_fast(cluster_map, w_ssd_map, b_ssd_map, q, pair):
    """
    Modify the maps in-place.
    The new cluster index is one greater than the max existing cluster index.
    @param cluster_map: maps a cluster index to a set of point indices
    @param w_ssd_map: maps a cluster index to a sum of squared distances
    @param b_ssd_map: maps a cluster index pair to a sum of squared distances
    @param q: a priority queue for near-neighbor cluster index pairs
    @param pair: the pair of cluster indices to merge
    """
    # define the cluster indices to be merged
    ia, ib = pair
    # define the index of the merged cluster
    neo = max(cluster_map) + 1
    # define the ssd keys that contain a constituent cluster index
    dead_ssd_keys = get_dead_ssd_keys(cluster_map, pair)
    # Add ssd to the new cluster to the between cluster sum of squares.
    # Also add these to the priority queue.
    for i in cluster_map:
        if i in pair:
            # do not compute distances to the constituent clusters
            continue
        ka = (i,ia) if i < ia else (ia,i)
        kb = (i,ib) if i < ib else (ib,i)
        b_ssd = b_ssd_map[ka] + b_ssd_map[kb]
        b_ssd_map[i,neo] = b_ssd
        heapq.heappush(q, (b_ssd, (i, neo)))
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

def get_pair_fast(cluster_map, q):
    """
    Pairs are popped from the queue until a pair whose indices
    are both current cluster indices is found.
    The pair at the front of the queue has smallest sum squared distance
    @param q: a priority queue
    @param cluster_map: a dict with cluster index keys
    """
    while True:
        d, (ia, ib) = heapq.heappop(q)
        if ia in cluster_map and ib in cluster_map:
            return ia, ib

def get_initial_queue(b_ssd_map):
    """
    @param points: a sequence of points in euclidean space
    @return: an initial priority queue
    """
    q = [(d, pair) for pair, d in b_ssd_map.items()]
    heapq.heapify(q)
    return q


class TestAgglom(unittest.TestCase):

    def test_placeholder(self):
        pass


if __name__ == '__main__':
    unittest.main()
