"""Cluster using k-means.
"""

from StringIO import StringIO
import os

import argparse
import numpy as np

from SnippetUtil import HandlingError
import Form
import Util
import Carbone
import kmeans

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
    wcss, labels = kmeans.lloyd_with_restarts(points, fs.k, nrestarts)
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
