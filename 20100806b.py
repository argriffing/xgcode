"""Find the number of agglomerated clusters using the gap statistic.

Apply hierarchical (agglomerative) clustering,
using squared error and average linkage.
This follows the protocol of Tibshirani et al. in example 4.1
of 'Estimating the number of clusters in a data set via the gap statistic'.
"""

from StringIO import StringIO
import os
import time
import random

import numpy as np

from SnippetUtil import HandlingError
import Form
import FormOut
import Util
import iterutils
import Carbone
import agglom
import kmeans
import tibshirani
import const
import RUtil

g_tags = ['pca:compute']

g_default = const.read('20100709a')

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('table', 'R table', g_default),
            Form.Sequence('axes', 'column labels of Euclidean axes',
                ('pc1', 'pc2', 'pc3')),
            Form.Integer('nsamples',
                'use this many samples for the null distribution', 10),
            Form.CheckGroup('options', 'more options', [
                Form.CheckItem('verbose',
                    'show index values', True)])]
    return form_objects

def get_form_out():
    """
    @return: the format of the output
    """
    return FormOut.Report('report')

def do_sampling(extents, npoints, B):
    """
    @param extents: the length of each axis of the hypercube
    @param npoints: sample this many points at a time
    @param B: use this many monte carlo samples
    @return: (k, expected log Wk, sk) triples
    """
    # Get the list of the number of clusters at each step,
    # and get the corresponding logarithms of wgss.
    wlog_arr = []
    for isample in range(B):
        pairs = get_nclusters_logw_pairs(extents, npoints)
        nclusters_list, wgss_list = zip(*pairs)
        wlog_arr.append(np.log(wgss_list))
    # Each row of this transpose
    # is a sample of wlogs for a given k clusters.
    W = np.array(wlog_arr).T
    # Get the expectations and the gap thresholds.
    expectations = []
    thresholds = []
    for wlogs in W:
        expectations.append(np.mean(wlogs))
        thresholds.append(tibshirani.get_simulation_correction(wlogs))
    if len(nclusters_list) != len(expectations):
        msg = 'expected as many expectations as cluster levels'
        raise ValueError(msg)
    if len(nclusters_list) != len(thresholds):
        msg = 'expected as many thresholds as cluster levels'
        raise ValueError(msg)
    # reverse all of the lists so that they are by increasing cluster size
    triples = list(reversed(zip(nclusters_list, expectations, thresholds)))
    # Return the nclusters_list, the expectations, and the thresholds.
    return zip(*triples)

def get_nclusters_logw_pairs(extents, npoints):
    """
    Get sample statistics by agglomeratively clustering random points.
    These sample statistics will be used for a null distribution.
    @param extents: the length of each axis of the hypercube
    @param npoints: sample this many points at a time
    @return: (nclusters, logw) pairs for a single sampling of points
    """
    # sample the points
    pointlist = []
    for i in range(npoints):
        p = [random.uniform(0, x) for x in extents]
        pointlist.append(p)
    points = np.array(pointlist)
    # do the clustering, recording the within group sum of squares
    nclusters_wgss_pairs = []
    allmeandist = kmeans.get_allmeandist(points)
    cluster_map = agglom.get_initial_cluster_map(points)
    b_ssd_map = agglom.get_initial_b_ssd_map(points)
    w_ssd_map = agglom.get_initial_w_ssd_map(points)
    q = agglom.get_initial_queue(b_ssd_map)
    while len(cluster_map) > 2:
        pair = agglom.get_pair_fast(cluster_map, q)
        agglom.merge_fast(cluster_map, w_ssd_map, b_ssd_map, q, pair)
        indices = cluster_map.keys()
        wgss = sum(w_ssd_map[i] / float(len(cluster_map[i])) for i in indices)
        nclusters_wgss_pairs.append((len(cluster_map), wgss))
    return nclusters_wgss_pairs

def get_response_content(fs):
    # read the table
    rtable = RUtil.RTable(fs.table.splitlines())
    header_row = rtable.headers
    data_rows = rtable.data
    Carbone.validate_headers(header_row)
    # get the numpy array of conformant points
    h_to_i = dict((h, i+1) for i, h in enumerate(header_row))
    axis_headers = fs.axes
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
            axis_list = Carbone.get_numeric_column(data_rows, index)
        except Carbone.NumericError:
            msg_a = 'expected the axis column %s ' % h
            msg_b = 'to be numeric'
            raise ValueError(msg_a + msg_b)
        axis_lists.append(axis_list)
    points = np.array(zip(*axis_lists))
    # do the clustering while computing the wgss at each merge
    cluster_counts = []
    wgss_values = []
    allmeandist = kmeans.get_allmeandist(points)
    cluster_map = agglom.get_initial_cluster_map(points)
    w_ssd_map = agglom.get_initial_w_ssd_map(points)
    b_ssd_map = agglom.get_initial_b_ssd_map(points)
    q = agglom.get_initial_queue(b_ssd_map)
    while len(cluster_map) > 2:
        # do an agglomeration step
        pair = agglom.get_pair_fast(cluster_map, q)
        agglom.merge_fast(cluster_map, w_ssd_map, b_ssd_map, q, pair)
        # compute the within group sum of squares
        indices = cluster_map.keys()
        wgss = sum(w_ssd_map[i] / float(len(cluster_map[i])) for i in indices)
        # compute the between group sum of squares
        bgss = allmeandist - wgss
        # append to the lists
        cluster_counts.append(len(cluster_map))
        wgss_values.append(wgss)
    # compute the log wgss values
    wlogs = np.log(wgss_values)
    # reverse the log values so that they are by increasing cluster size
    wlogs = list(reversed(wlogs))
    # sample from the null distribution
    extents = np.max(points, axis=0) - np.min(points, axis=0)
    nclusters_list, expectations, thresholds = do_sampling(
            extents, len(points), fs.nsamples)
    # get the gaps
    gaps = np.array(expectations) - wlogs
    # Get the best cluster count according to the gap statistic.
    best_i = None
    criteria = []
    for i, ip1 in iterutils.pairwise(range(len(nclusters_list))):
        k, kp1 = nclusters_list[i], nclusters_list[ip1]
        criterion = gaps[i] - gaps[ip1] + thresholds[ip1]
        criteria.append(criterion)
        if criterion > 0:
            if best_i is None:
                best_i = i
    best_k = nclusters_list[best_i]
    # create the response
    out = StringIO()
    print >> out, 'best cluster count: k = %d' % best_k
    if fs.verbose:
        print >> out
        print >> out, '(k, expected, observed, gap, threshold, criterion):'
        n = len(nclusters_list)
        for i, k in enumerate(nclusters_list):
            row = [k, expectations[i], wlogs[i], gaps[i], thresholds[i]]
            if i < n-1:
                row += [criteria[i]]
            else:
                row += ['-']
            print >> out, '\t'.join(str(x) for x in row)
    # return the response
    return out.getvalue()
