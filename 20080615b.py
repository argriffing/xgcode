"""Convert a path resistance matrix to a laplacian matrix.
"""

import numpy

from SnippetUtil import HandlingError
import SnippetUtil
import MatrixUtil
import Euclid
import Form

def get_form():
    """
    @return: the body of a form
    """
    D = numpy.array([
            [0, 2, 2],
            [2, 0, 2],
            [2, 2, 0]])
    return [Form.Matrix('matrix', 'path resistance matrix', D, MatrixUtil.assert_predistance)]

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    D = fs.matrix
    L = Euclid.edm_to_laplacian(D)
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, MatrixUtil.m_to_string(L)
