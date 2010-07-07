"""
Use the standard k-means clustering algorithm.
"""

import random

import numpy as np

import iterutils
import MatrixUtil


def get_point_center_sqdists(points, centers):
    """
    Inputs and outputs are numpy arrays.
    @param points: the points to be reclustered
    @param centers: cluster centers
    @return: for each point, the squared distance to each center
    """
    MatrixUtil.assert_2d(points)
    MatrixUtil.assert_2d(centers)
    npoints = len(points)
    ncenters = len(centers)
    # get the dot products of points with themselves
    #pself = np.sum(points*points, axis=1)
    pself = np.array([np.dot(p, p) for p in points])
    # get the dot products of centers with themselves
    #cself = np.sum(centers*centers, axis=1)
    cself = np.array([np.dot(c, c) for c in centers])
    # get the matrix product of points and centers
    prod = np.dot(points, centers.T)
    # get the matrix of squared distances
    sqdists = (
        np.outer(pself, np.ones(ncenters)) +
        np.outer(np.ones(npoints), cself) -
        2*prod)
    return sqdists

def get_centers(points, labels):
    """
    Inputs and outputs are numpy arrays.
    @param points: euclidean points
    @param labels: conformant cluster indices
    """
    MatrixUtil.assert_2d(points)
    MatrixUtil.assert_1d(labels)
    if len(points) != len(labels):
        raise ValueError('array incompatibility')
    ncoords = len(points[0])
    nclusters = max(labels) + 1
    sums = [np.zeros(ncoords) for i in range(nclusters)]
    counts = [0]*nclusters
    for point, label in zip(points, labels):
        sums[label] += point
        counts[label] += 1
    M = np.array([s/c for s, c in zip(sums, counts)])
    return M

def get_labels(sqdists):
    """
    Inputs and outputs are numpy arrays.
    Account for the fact that sometimes a cluster will go away.
    That is, if no point is in the voronoi region of a centroid,
    then in the next iteration this cluster should disappear.
    @param sqdists: for each point, the squared distance to each center
    @return: for each point, the label of the nearest cluster
    """
    labels = np.argmin(sqdists, axis=1)
    new_to_old = list(iterutils.unique_everseen(labels))
    old_to_new = dict((old, new) for new, old in enumerate(new_to_old))
    return np.array([old_to_new[old] for old in labels])

def get_wcss(sqdists, labels):
    """
    Get the within-cluster sum of squares.
    @param sqdists: for each point, the squared distance to each center
    @param labels: cluster labels
    @return: within-cluster sum of squares
    """
    MatrixUtil.assert_2d(sqdists)
    MatrixUtil.assert_1d(labels)
    if len(sqdists) != len(labels):
        raise ValueError('array incompatibility')
    return sum(row[label] for row, label in zip(sqdists, labels))

def get_random_labels(npoints, nclusters):
    """
    Get random labels with each label appearing at least once.
    """
    labels = np.random.randint(0, nclusters, npoints)
    indices = random.sample(range(npoints), nclusters)
    for index, label in zip(indices, range(nclusters)):
        labels[index] = label
    return labels

def lloyd(points, labels):
    """
    This is the standard algorithm for kmeans clustering.
    @param points: points in euclidean space
    @param labels: initial cluster labels
    @return: within cluster sum of squares, and labels
    """
    while True:
        centers = get_centers(points, labels)
        sqdists = get_point_center_sqdists(points, centers)
        next_labels = get_labels(sqdists)
        if np.array_equal(next_labels, labels):
            wcss = get_wcss(sqdists, labels)
            return wcss, labels
        labels = next_labels

def lloyd_with_restarts(points, nclusters, nrestarts):
    """
    This is the standard algorithm for kmeans clustering with restarts.
    @param points: points in euclidean space
    @param nclusters: the number of clusters
    @param nrestarts: the number of random restarts
    @return: within cluster sum of squares, labels
    """
    npoints = len(points)
    best_wcss = None
    best_labels = None
    for i in range(nrestarts):
        labels = get_random_labels(npoints, nclusters)
        wcss, labels = lloyd(points, labels)
        if (best_wcss is None) or (wcss < best_wcss):
            best_wcss = wcss
            best_labels = labels
    return best_wcss, best_labels
