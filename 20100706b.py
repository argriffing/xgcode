"""Cluster using k-means.

Input is an R table.
Output is an R table with an extra column.
"""

from StringIO import StringIO
import os
import time

import argparse
import numpy as np

from SnippetUtil import HandlingError
import Form
import FormOut
import Util
import Carbone
import kmeans
import combobreaker
import const
import moretypes

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
            Form.Integer('k', 'maximum number of clusters',
                2, low=2),
            Form.SingleLine('annotation', 'header of added column',
                'cluster'),
            kmeans.InitStrategy(),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.RTable('out')

def get_rtable_info(rtable, cluster_header, axis_headers):
    """
    @param rtable: a Carbone.RTable object
    @param cluster_header: header of the new column to add
    @param axis_headers: a tuple of column headers
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
    nrestarts = 10
    init_strategy = kmeans.InitStrategy().string_to_function(fs.kmeans_init)
    wcss, labels = kmeans.lloyd_with_restarts(
            points, fs.k, nrestarts, init_strategy)
    # get the response
    lines = ['\t'.join(header_row + [fs.annotation])]
    for i, (label, data_row) in enumerate(zip(labels, data_rows)):
        row = data_row + [str(label)]
        lines.append('\t'.join(row))
    # return teh response
    return '\n'.join(lines) + '\n'


class GlobalState(object):
    def __init__(self, rtable, points, annotation, nclusters, init_strategy):
        """
        @param rtable: a Carbone.RTable object
        @param points: points as rows in a numpy array
        @param annotation: the name of the new column header
        @param nclusters: target this many clusters
        @param init_strategy: use this strategy to guess initial centers
        """
        self.rtable = rtable
        self.points = points
        self.annotation = annotation
        self.nclusters = nclusters
        self.init_strategy = init_strategy

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
        centers = np.array(list(self.init_strategy(
            self.gs.points, self.gs.nclusters)))
        sqdists = kmeans.get_point_center_sqdists(self.gs.points, centers)
        labels = kmeans.get_labels(sqdists)
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


def gen_states(gs):
    """
    @param gs: state that is unchanged between iterations
    """
    state = ClusterState(gs)
    while True:
        state = state.get_next()
        yield state


def main(args):
    """
    @param args: from argparse
    """
    # get some state that will not change between k-means restarts
    with open(args.table_filename) as fin:
        rtable = Carbone.RTable(fin)
    points = get_rtable_info(rtable, args.annotation, args.axes)
    init_strategy = kmeans.InitStrategy().string_to_function(args.kmeans_init)
    gs = GlobalState(rtable, points, args.annotation, args.k, init_strategy)
    # go until iteration is stopped for some reason
    run_info = combobreaker.combo_breaker(
            gen_states(gs), args.nseconds, args.nrestarts)
    print run_info.get_response()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--table_filename', required=True,
            help='R table filename')
    parser.add_argument('--axes', required=True,
            help='column labels of Euclidean axes')
    parser.add_argument('--k', type=moretypes.int_ge_2, required=True,
            help='target number of clusters')
    parser.add_argument('--annotation', default='cluster',
            help='header of added column')
    parser.add_argument('--nrestarts', type=moretypes.whole_number,
            help='restart the k-means iterative refinement this many times')
    parser.add_argument('--nseconds', type=moretypes.positive_float,
            help='run for this many seconds')
    kmeans.InitStrategy().add_argument(parser)
    main(parser.parse_args())
