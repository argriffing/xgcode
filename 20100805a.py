"""Cluster using agglomerative clustering.

Input is an R table.
Output is an R table with an extra column.
Apply hierarchical (agglomerative) clustering,
using squared error and average linkage.
This follows the protocol of Tibshirani et al. in example 4.1
of 'Estimating the number of clusters in a data set via the gap statistic'.
"""

import numpy as np

from SnippetUtil import HandlingError
import Form
import FormOut
import Util
import Carbone
import agglom
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
            Form.Integer('k', 'maximum number of clusters',
                2, low=2),
            Form.SingleLine('annotation', 'header of added column',
                'cluster'),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.RTable('out')

def get_rtable_info(rtable, cluster_header, axes_string):
    """
    @param rtable: a Carbone.RTable object
    @param cluster_header: header of the new column to add
    @param axes_string: something like "pc1 pc2 pc2" with column headers
    @return: points as rows in a numpy array
    """
    header_row = rtable.headers
    data_rows = rtable.data
    # do header validation
    Carbone.validate_headers(header_row)
    if not Carbone.is_valid_header(cluster_header):
        raise ValueError('invalid column header: %s' % cluster_header)
    if cluster_header in header_row:
        msg = 'the column header %s is already in the table' % cluster_header
        raise ValueError(msg)
    # get the numpy array of conformant points
    h_to_i = dict((h, i+1) for i, h in enumerate(header_row))
    axis_headers = axes_string.split()
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
    return points

def get_response_content(fs):
    rtable = Carbone.RTable(fs.table.splitlines())
    header_row = rtable.headers
    data_rows = rtable.data
    points = get_rtable_info(rtable, fs.annotation, fs.axes)
    # do the clustering
    cluster_map = agglom.get_initial_cluster_map(points)
    w_ssd_map = agglom.get_initial_w_ssd_map(points)
    b_ssd_map = agglom.get_initial_b_ssd_map(points)
    q = agglom.get_initial_queue(b_ssd_map)
    while len(cluster_map) > fs.k:
        pair = agglom.get_pair_fast(cluster_map, q)
        agglom.merge_fast(cluster_map, w_ssd_map, b_ssd_map, q, pair)
    # create the map from a point index to a cluster index
    point_to_cluster = {}
    for cluster_index, point_indices in cluster_map.items():
        for point_index in point_indices:
            point_to_cluster[point_index] = cluster_index
    # define the raw labels which may be big numbers
    raw_labels = [point_to_cluster[i] for i, p in enumerate(points)]
    # rename the labels with small numbers
    raw_to_label = dict((b, a) for  a, b in enumerate(sorted(set(raw_labels))))
    labels = [raw_to_label[raw] for raw in raw_labels]
    # get the response
    lines = ['\t'.join(header_row + [fs.annotation])]
    for i, (label, data_row) in enumerate(zip(labels, data_rows)):
        row = data_row + [str(label)]
        lines.append('\t'.join(row))
    # return the response
    return '\n'.join(lines) + '\n'
