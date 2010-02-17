"""Create some example graph matrices and layout data.
"""


from StringIO import StringIO
import os
import math

import numpy as np

from SnippetUtil import HandlingError
import SnippetUtil
import Form
import GPS
import Euclid

def get_locations():
    sqrt2 = math.sqrt(2.0)
    points = [
            (-1-sqrt2, 1),
            (-1-sqrt2, -1),
            (1+sqrt2, 1),
            (4+sqrt2, -4),
            (-sqrt2, 0),
            (sqrt2, 0)]
    return points

def get_edges():
    edges = [
            (0, 4),
            (1, 4),
            (2, 5),
            (3, 5),
            (4, 5)]
    return edges

def get_form():
    """
    @return: a list of form objects
    """
    # define the list of form objects
    form_objects = []
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    locations = get_locations()
    np_locs = [np.array(p) for p in locations]
    edges = get_edges()
    npoints = len(locations)
    # start writing the response
    np.set_printoptions(linewidth=200)
    out = StringIO()
    # print the layout data
    print >> out, 'POINTS'
    for i, (x, y) in enumerate(locations):
        print >> out, i, x, y
    print >> out, 'EDGES'
    for i, j in edges:
        print >> out, i, j
    print >> out
    # show the unweighted adjacency matrix
    UA = np.zeros((npoints, npoints))
    for i, j in edges:
        UA[i, j] = 1
        UA[j, i] = 1
    print >> out, 'unweighted adjacency matrix:'
    print >> out, UA
    print >> out
    # show the unweighted laplacian matrix
    UL = Euclid.adjacency_to_laplacian(UA)
    print >> out, 'unweighted laplacian matrix:'
    print >> out, UL
    print >> out
    # show the weighted adjacency matrix
    WA = np.zeros((npoints, npoints))
    for i, j in edges:
        d = np.linalg.norm(np_locs[i] - np_locs[j]) / math.sqrt(2.0)
        w = 1.0 / d
        WA[i, j] = w
        WA[j, i] = w
    print >> out, 'weighted adjacency matrix:'
    print >> out, WA
    print >> out
    # show the weighted laplacian matrix
    WL = Euclid.adjacency_to_laplacian(WA)
    print >> out, 'weighted laplacian matrix:'
    print >> out, WL
    print >> out
    # write the response
    return [('Content-Type', 'text/plain')], out.getvalue().strip()
