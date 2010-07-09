"""Cluster using k-means.
"""

from StringIO import StringIO
import os
import time

import argparse
import numpy as np

from SnippetUtil import HandlingError
import Form
import Util
import Carbone
import kmeans
import combobreaker

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
            axis_list = get_numeric_column(data_rows, index)
        except NumericError:
            msg_a = 'expected the axis column %s ' % h
            msg_b = 'to be numeric'
            raise ValueError(msg_a + msg_b)
        axis_lists.append(axis_list)
    points = np.array(zip(*axis_lists))
    return points

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    rtable = Carbone.RTable(fs.table.splitlines())
    header_row = rtable.headers
    data_rows = rtable.data
    points = get_rtable_info(rtable, fs.annotation, fs.axes)
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


class GlobalState(object):
    def __init__(self, rtable, points, annotation, nclusters):
        """
        @param rtable: a Carbone.RTable object
        @param points: points as rows in a numpy array
        @param annotation: the name of the new column header
        @param nclusters: target this many clusters
        """
        self.rtable = rtable
        self.points = points
        self.annotation = annotation
        self.nclusters = nclusters

class ClusterState(object):
    # TODO this should possibly be a real generator

    def __init__(self, global_state):
        """
        @param global_state: things that do not change between iterations
        """
        self.gs = global_state
        self.best_wcss = None
        self.best_labels = None

    def get_response(self):
        """
        @return: a string
        """
        if self.best_wcss is None:
            return 'no results'
        else:
            lines = ['\t'.join(self.gs.rtable.headers + [self.gs.annotation])]
            label_row_pairs = zip(self.best_labels, self.gs.rtable.data)
            for i, (label, data_row) in enumerate(label_row_pairs):
                row = data_row + [str(label)]
                lines.append('\t'.join(row))
            return '\n'.join(lines)

    def mutate(self):
        """
        Do an iteration of the Lloyd algorithm.
        """
        npoints = len(gs.points)
        labels = kmeans.get_random_labels(npoints, self.gs.nclusters)
        wcss, labels = kmeans.lloyd(self.gs.points, labels)
        if (self.best_wcss is None) or (wcss < self.best_wcss):
            self.best_wcss = wcss
            self.best_labels = labels

    def get_next(self):
        """
        This call should return a new object without self-mutation.
        """
        next_state = ClusterState(self.gs)
        next_state.best_wcss = self.best_wcss
        next_state.best_labels = self.best_labels
        next_state.mutate()
        return next_state


def get_next_state(state):
    return state.get_next()

def main(args):
    """
    @param args: argparse'd
    """
    # get some state that will not change between k-means restarts
    with open(args.table_filename) as fin:
        rtable = Carbone.RTable(fin)
    points = get_rtable_info(rtable, args.annotation, args.axes)
    gs = GlobalState(rtable, points, args.annotation, args.k)
    # go until iteration is stopped for some reason
    run_info = combobreaker.combo_breaker(
            ClusterState(gs), get_next_state, args.nseconds, args.nrestarts)
    print run_info.get_response()

if __name__ == '__main__':
    parser = arparse.ArgumentParser(description=__doc__)
    parser.add_argument('--table_filename', required=True,
            help='R table filename')
    parser.add_argument('--axes', required=True,
            help='column labels of Euclidean axes')
    parser.add_argument('--k', type=int_ge_2, required=True,
            help='target number of clusters')
    parser.add_argument('--annotation', default='cluster',
            help='header of added column')
    parser.add_argument('--nrestarts', type=whole_number,
            help='restart the k-means iterative refinement this many times')
    parser.add_argument('--nseconds', type=positive_float,
            help='run for this many seconds')
    main(parser.parse_args())

def whole_number(x):
    try:
        x = int(x)
    except ValueError:
        raise TypeError
    if x < 1:
        raise TypeError
    return x

def positive_float(x):
    try:
        x = float(x)
    except ValueError:
        raise TypeError
    if x < 0:
        raise TypeError
    return x

def int_ge_2(x):
    x = whole_number(x)
    if x < 2:
        raise TypeError
    return x

