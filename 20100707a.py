"""Find the number of clusters for k-means using the Calinski-Harabasz index.
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
            Form.CheckGroup('options', 'more options', [
                Form.CheckItem('verbose',
                    'show calinski index values', True)])]
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
            axis_list = get_numeric_column(data_rows, index)
        except NumericError:
            msg_a = 'expected the axis column %s ' % h
            msg_b = 'to be numeric'
            raise ValueError(msg_a + msg_b)
        axis_lists.append(axis_list)
    points = np.array(zip(*axis_lists))
    # precompute some stuff
    allmeandist = get_allmeandist(points)
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
        wgss, labels = kmeans.lloyd_with_restarts(points, k, nrestarts)
        bgss = allmeandist - wgss
        calinski = get_calinski_index(bgss, wgss, k, n)
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
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue()

def get_allmeandist(points):
    """
    @param points: each row is a point and each column is a coordinate
    @return: a single number
    """
    return np.sum((points - np.mean(points, axis=0))**2)

def get_calinski_index(bgss, wgss, k, n):
    """
    Choose the number of clusters that gives the greatest calinski index.
    @param bgss: between groups sum of squares
    @param wgss: within groups sum of squares
    @param k: number of clusters
    @param n: number of points
    @return: a floating point number
    """
    if not (1 < k < n):
        msg_a = 'the calinski index '
        msg_b = 'is defined for integers k and n such that 1 < k < n'
        raise ValueError(msg_a + msg_b)
    numerator = bgss / float(k - 1)
    denominator = wgss / float(n - k)
    if not denominator:
        return float('inf')
    return numerator / denominator

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
