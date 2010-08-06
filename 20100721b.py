"""Cluster using the scipy k-means implementation.

Input is an R table.
Output is verbose commentary and debugging information.
"""

from StringIO import StringIO
import os
import time

import argparse
import numpy as np
from scipy import cluster

from SnippetUtil import HandlingError
import Form
import FormOut
import Util
import Carbone
import kmeans
import combobreaker
import const

g_tags = ['pca:misc']

g_default_lines = [
        'p1 p2',
        '1 -20 -20',
        '2 -30 -20',
        '3 20 20',
        '4 20 26']

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('table', 'R table', '\n'.join(g_default_lines)),
            Form.SingleLine('axes', 'column labels of Euclidean axes',
                ' '.join(('p1', 'p2'))),
            Form.Integer('k', 'maximum number of clusters',
                2, low=2)]
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
    # define constants
    nrestarts = 10
    # read the input
    rtable = Carbone.RTable(fs.table.splitlines())
    header_row = rtable.headers
    data_rows = rtable.data
    points = get_rtable_info(rtable, 'annotation', fs.axes)
    # do the clustering
    codebook, distortion = cluster.vq.kmeans(
            points, fs.k, iter=nrestarts, thresh=1e-9)
    sqdists = kmeans.get_point_center_sqdists(points, codebook)
    labels = kmeans.get_labels_without_cluster_removal(sqdists)
    wgss = kmeans.get_wcss(sqdists, labels)
    norms = [np.linalg.norm(p-codebook[g]) for p, g in zip(points, labels)]
    redistortion = np.mean(norms)
    # create the response
    out = StringIO()
    print >> out, 'scipy distortion:', distortion
    print >> out, 'recomputed distortion:', redistortion
    print >> out, 'wgss:', wgss
    # return the response
    return out.getvalue()
