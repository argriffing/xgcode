"""Normalize a rate matrix to have one expected transition per unit time.
"""

from StringIO import StringIO

import numpy

from SnippetUtil import HandlingError
import MatrixUtil
import RateMatrix
import Form

def get_form():
    """
    @return: the body of a form
    """
    # define the default matrix string
    R = numpy.array([
            [-3, 1, 1, 1],
            [1, -3, 1, 1],
            [1, 1, -3, 1],
            [1, 1, 1, -3]])
    # define the form objects
    return [Form.Matrix('matrix', 'rate matrix', R, MatrixUtil.assert_rate_matrix)]

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # read the matrix from the form data
    R = fs.matrix
    n = len(R)
    # convert the row major rate matrix to a rate matrix object
    arbitrary_states = [str(x) for x in range(n)]
    rate_matrix_object = RateMatrix.RateMatrix(R.tolist(), arbitrary_states)
    rate_matrix_object.normalize()
    normalized_row_major = rate_matrix_object.get_row_major_rate_matrix()
    rate_matrix_string = MatrixUtil.m_to_string(normalized_row_major)
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, rate_matrix_string
