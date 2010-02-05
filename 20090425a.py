"""Do varimax rotation.
"""

from StringIO import StringIO

import numpy

from SnippetUtil import HandlingError
import Util
import MatrixUtil
import Form

def get_form():
    """
    @return: the body of a form
    """
    return [Form.Matrix('matrix', 'matrix', MatrixUtil.g_example_loading_matrix)]

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    rotated_matrix = numpy.dot(fs.matrix, MatrixUtil.varimax(fs.matrix))
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, MatrixUtil.m_to_string(rotated_matrix)
