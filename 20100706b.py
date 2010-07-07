"""Cluster using k-means.
"""

from StringIO import StringIO
import os
import random

import argparse
import numpy as np

from SnippetUtil import HandlingError
import Form
import Util
import Carbone
import iterutils

g_default_rows = [
        ['otu', 'species', 'location', 'temperature', 'precipitation',
            'pc1', 'pc2', 'pc3'],
        [1, 'IC100', 'Ap', 'GA', 15.0, 600.0,
            -2.8053476259, 0.556532380058, -6.17891756957],
        [2, 'IC101', 'Ap', 'GA', 15.0, 600.0,
            -2.8053476259, 0.556532380058, -6.17891756956],
        [3, 'IC102', 'Ap', 'GA', 15.0, 600.0,
            -2.80534762591, 0.556532380059, -6.17891756957],
        [455, 'IC577', 'Ac', 'AR', 25.0, 400.0,
            -13.7544736082, -7.16259232881, 7.0902951321],
        [456, 'IC580', 'Ac', 'AR', 25.0, 400.0,
            3.56768959361, 0.385873934264, 1.23981735331],
        [457, 'IC591', 'Ac', 'AR', 25.0, 400.0,
            -11.6455270418, -5.710582374, 5.60835091179]]

g_default_lines = ['\t'.join(str(x) for x in row) for row in g_default_rows]
g_default_string = '\n'.join(g_default_lines)

def get_point_center_sqdists(points, centers):
    """
    Inputs and outputs are numpy arrays.
    @param points: the points to be reclustered
    @param centers: cluster centers
    @return: for each point, the squared distance to each center
    """
    if len(points.shape) != 2:
        raise ValueError('expected a point matrix')
    if len(centers.shape) != 2:
        raise ValueError('expected a matrix of cluster centers')
    npoints = len(points)
    ncenters = len(centers)
    # get the dot products of points with themselves
    pself = np.array([np.dot(p,p) for p in points])
    # get the dot products of centers with themselves
    cself = np.array([np.dot(c,c) for c in centers])
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
    ncoords = len(points[0])
    nclusters = max(labels) + 1
    sums = [np.zeros(ncoords) for i in range(nclusters)]
    counts = [0]*nclusters
    for point, label in zip(points, labels):
        sums[label] += point
        counts[label] += 1
    M = np.array([s/c for s, c in zip(sums, counts)])
    return M

def get_normalized_labels(labels):
    """
    Account for the fact that sometimes a cluster will go away.
    That is, if no point is in the voronoi region of a centroid,
    then in the next iteration this cluster should disappear.
    @param labels: cluster labels
    @return: cluster labels
    """
    new_to_old = list(iterutils.unique_everseen(labels))
    old_to_new = dict((old, new) for new, old in enumerate(new_to_old))
    return np.array([old_to_new[old] for old in labels])

def get_labels(sqdists):
    """
    Inputs and outputs are numpy arrays.
    @param sqdists: for each point, the squared distance to each center
    @return: for each point, the label of the nearest cluster
    """
    return get_normalized_labels(np.argmin(sqdists, axis=1))

def get_wcss(sqdists, labels):
    """
    Get the within-cluster sum of squares.
    @param sqdists: for each point, the squared distance to each center
    @param labels: cluster labels
    @return: within-cluster sum of squares
    """
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
    @return: labels
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
    return best_labels


def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('table', 'R table', g_default_string),
            Form.SingleLine('axes', 'column labels of Euclidean axes',
                ' '.join(('pc1', 'pc2', 'pc3'))),
            Form.Integer('k', 'maximum number of clusters',
                2, low=2),
            Form.SingleLine('annotation', 'header of added column',
                'cluster'),
            Form.ContentDisposition()]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # read the table
    rtable = Carbone.RTable(fs.table.splitlines())
    header_row = rtable.headers
    data_rows = rtable.data
    # do header validation
    Carbone.validate_headers(header_row)
    if not Carbone.is_valid_header(fs.annotation):
        raise ValueError('invalid column header: %s' % fs.annotation)
    if fs.annotation in header_row:
        msg = 'the column header %s is already in the table' % fs.annotation
        raise ValueError(msg)
    # get the numpy array of conformant points
    h_to_i = dict((h, i+1) for i, h in enumerate(header_row))
    axis_headers = fs.axes.split()
    if not axis_headers:
        raise ValueError('no Euclidean axes were provided')
    axis_set = set(axis_headers)
    header_set = set(header_row)
    bad_axes = axis_set - header_set
    if bad_axes:
        raise ValueError('invalid axes: ' + ', '.join(bad_axes))
    axis_lists = []
    for h in axis_headers:
        index = h_to_i[h]
        try:
            axis_list = get_numeric_column(data_rows, index)
        except NumericError:
            msg_a = 'expected the axis column %s ' % h
            msg_b = 'to be numeric'
            raise ValueError(msg_a + msg_b)
        axis_lists.append(axis_list)
    points = np.array(zip(*axis_lists))
    # do the clustering
    nrestarts = 10
    labels = lloyd_with_restarts(points, fs.k, nrestarts)
    # get the response
    lines = ['\t'.join(header_row + [fs.annotation])]
    for i, (label, data_row) in enumerate(zip(labels, data_rows)):
        row = data_row + [str(label)]
        lines.append('\t'.join(row))
    content = '\n'.join(lines) + '\n'
    # return the response
    disposition = '%s; filename=%s' % (fs.contentdisposition, 'out.table')
    response_headers = [
            ('Content-Type', 'text/plain'),
            ('Content-Disposition', disposition)]
    return response_headers, content

class NumericError(Exception): pass

def get_numeric_column(data, index):
    """
    @param data: row major list of lists of numbers as strings
    @param index: column index
    @return: list of floats
    """
    strings = zip(*data)[index]
    try:
        floats = [float(x) for x in strings]
    except ValueError, v:
        raise NumericError
    return floats
