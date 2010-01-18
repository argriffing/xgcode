"""Given a weighted adjacency matrix, calculate the path resistance matrix.

Here the weighted adjacency matrix is the conductance matrix.
Each element is the reciprocal of the value of a resistor directly connecting the nodes.
"""

import StringIO

import numpy

from SnippetUtil import HandlingError
import MatrixUtil
import Euclid
import Form

def get_form():
    """
    @return: the body of a form
    """
    A = numpy.array([
        [0, 2, 2, 0, 0, 0, 0],
        [2, 0, 2, 0, 0, 0, 0],
        [2, 2, 0, 3, 0, 0, 0],
        [0, 0, 3, 0, 2, 2, 0],
        [0, 0, 0, 2, 0, 2, 1],
        [0, 0, 0, 2, 2, 0, 1],
        [0, 0, 0, 0, 1, 1, 0]])
    return [Form.Matrix('matrix', 'weighted adjacency matrix', A, MatrixUtil.assert_weighted_adjacency)]

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    A = fs.matrix
    L = Euclid.adjacency_to_laplacian(A)
    D = Euclid.laplacian_to_edm(L)
    # start to prepare the reponse
    out = StringIO.StringIO()
    print >> out, MatrixUtil.m_to_string(D)
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
