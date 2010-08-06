"""Find the number of agglomerated clusters using the Calinski-Harabasz index.

Apply hierarchical (agglomerative) clustering,
using squared error and average linkage.
This follows the protocol of Tibshirani et al. in example 4.1
of 'Estimating the number of clusters in a data set via the gap statistic'.
"""

from StringIO import StringIO
import os
import time

import numpy as np

from SnippetUtil import HandlingError
import Form
import FormOut
import Util
import Carbone
import agglom
import kmeans
import const

g_tags = ['pca:compute']

g_default = const.read('20100709a')

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('table', 'R table', g_default),
            Form.SingleLine('axes', 'column labels of Euclidean axes',
                ' '.join(('pc1', 'pc2', 'pc3'))),
            Form.CheckGroup('options', 'more options', [
                Form.CheckItem('verbose',
                    'show calinski index values', True)])]
    return form_objects

def get_form_out():
    """
    @return: the format of the output
    """
    return FormOut.Report('report')

def get_response_content(fs):
    # read the table
    rtable = Carbone.RTable(fs.table.splitlines())
    header_row = rtable.headers
    data_rows = rtable.data
    Carbone.validate_headers(header_row)
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
            axis_list = Carbone.get_numeric_column(data_rows, index)
        except Carbone.NumericError:
            msg_a = 'expected the axis column %s ' % h
            msg_b = 'to be numeric'
            raise ValueError(msg_a + msg_b)
        axis_lists.append(axis_list)
    points = np.array(zip(*axis_lists))
    # do the clustering while computing the calinski index at each merge
    cluster_counts = []
    wgss_values = []
    neg_calinskis = []
    allmeandist = kmeans.get_allmeandist(points)
    cluster_map = agglom.get_initial_cluster_map(points)
    w_ssd_map = agglom.get_initial_w_ssd_map(points)
    b_ssd_map = agglom.get_initial_b_ssd_map(points)
    while len(cluster_map) > 2:
        # do an agglomeration step
        pair = agglom.get_pair(cluster_map, b_ssd_map)
        agglom.merge(cluster_map, w_ssd_map, b_ssd_map, pair)
        # compute the within group sum of squares
        indices = cluster_map.keys()
        wgss = sum(w_ssd_map[i] / float(len(cluster_map[i])) for i in indices)
        # compute the between group sum of squares
        bgss = allmeandist - wgss
        # get the calinksi index
        n = len(points)
        k = len(cluster_map)
        numerator = bgss / float(k - 1)
        denominator = wgss / float(n - k)
        calinski = numerator / denominator
        # append to the lists
        cluster_counts.append(k)
        wgss_values.append(wgss)
        neg_calinskis.append(-calinski)
    # Get the best cluster count according to the calinski index.
    # Do this trickery with negs so that it breaks ties
    # using the smallest number of clusters.
    neg_calinksi, best_k = min(zip(neg_calinskis, cluster_counts))
    # create the response
    out = StringIO()
    print >> out, 'best cluster count: k = %d' % best_k
    if fs.verbose:
        print >> out
        print >> out, '(k, wgss, calinski):'
        for k, wgss, neg_calinski in zip(
                cluster_counts, wgss_values, neg_calinskis):
            row = (k, wgss, -neg_calinski)
            print >> out, '\t'.join(str(x) for x in row)
    # return the response
    return out.getvalue()
