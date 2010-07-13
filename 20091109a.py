"""Create a matrix representation of a path graph with equal weights.
"""

from StringIO import StringIO

import numpy as np

from SnippetUtil import HandlingError
import MatrixUtil
import Euclid
import iterutils
import Form
import FormOut

def get_form():
    """
    @return: a list of form objects
    """
    form_objects = [
            Form.Float('edge_affinity', 'edge affinities',
                1.0, low_exclusive=0),
            Form.Integer('nvertices', 'number of vertices',
                5, low=0, high=100),
            Form.RadioGroup('format', 'output options', [
                Form.RadioItem('adjacency', 'adjacency matrix'),
                Form.RadioItem('laplacian', 'laplacian matrix', True)])]
    return form_objects

def get_form_out():
    return FormOut.Matrix()

def create_adjacency_matrix(affinity, nvertices):
    """
    @param affinity: affinity between adjacent vertices
    @param nvertices: the number of vertices in the graph
    @return: a numpy matrix
    """
    A = np.zeros((nvertices, nvertices))
    for i, j in iterutils.pairwise(range(nvertices)):
        A[i,j] = affinity
        A[j,i] = affinity
    return A

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # create the matrix
    M = create_adjacency_matrix(fs.edge_affinity, fs.nvertices)
    if fs.laplacian:
        M = Euclid.adjacency_to_laplacian(M)
    # create the response
    out = StringIO()
    print >> out, MatrixUtil.m_to_string(M)
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
