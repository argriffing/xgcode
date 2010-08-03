"""Find the number of clusters for k-means using the Calinski-Harabasz index.
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
import const

g_tags = ['carbone_lab']

g_default = const.read('20100709a')

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('table', 'R table', g_default),
            Form.SingleLine('axes', 'column labels of Euclidean axes',
                ' '.join(('pc1', 'pc2', 'pc3'))),
            kmeans.InitStrategy(),
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
    # get the initialization strategy
    init_strategy = kmeans.InitStrategy().string_to_function(fs.kmeans_init)
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
    # precompute some stuff
    allmeandist = kmeans.get_allmeandist(points)
    nrestarts = 10
    nseconds = 2
    tm = time.time()
    n = len(points)
    wgss_list = []
    # neg because both items in the pair are used for sorting
    neg_calinski_k_pairs = []
    # look for the best calinski index in a small amount of time
    k = 2
    while True:
        wgss, labels = kmeans.lloyd_with_restarts(
                points, k, nrestarts, init_strategy)
        bgss = allmeandist - wgss
        calinski = kmeans.get_calinski_index(bgss, wgss, k, n)
        k_unique = len(set(labels))
        neg_calinski_k_pairs.append((-calinski, k_unique))
        wgss_list.append(wgss)
        if time.time() - tm > nseconds:
            break
        if k == n-1:
            break
        k += 1
    max_k = k
    best_neg_calinski, best_k = min(neg_calinski_k_pairs)
    best_calinski = -best_neg_calinski
    # create the response
    out = StringIO()
    print >> out, 'best cluster count: k = %d' % best_k
    print >> out, 'searched 2 <= k <= %d clusters' % max_k
    print >> out, '%.2f seconds' % (time.time() - tm)
    if fs.verbose:
        print >> out
        print >> out, '(k_unique, wgss, calinski):'
        for wgss, neg_calinski_k_pair in zip(wgss_list, neg_calinski_k_pairs):
            neg_calinski, k_unique = neg_calinski_k_pair
            calinski = -neg_calinski
            row = [k_unique, wgss, calinski]
            print >> out, '\t'.join(str(x) for x in row)
    # return the response
    return out.getvalue()

def main(args):
    #TODO something
    pass

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--table_filename', required=True,
            help='R table filename')
    parser.add_argument('--axes', required=True,
            help='column labels of Euclidean axes')
    parser.add_argument('--k_limit', type=int_ge_2,
            help='limit the search to at most this many clusters')
    parser.add_argument('--nseconds', type=positive_float,
            help='run for this many seconds')
    main(parser.parse_args())
