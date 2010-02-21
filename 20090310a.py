"""Given a distance matrix, get a covariance-like neighbor joining matrix.

Given a distance matrix,
get the covariance matrix implied by neighbor joining.
"""

from StringIO import StringIO

import numpy as np

from SnippetUtil import HandlingError
import MatrixUtil
import Form

def get_form():
    """
    @return: the body of a form
    """
    D = np.array([
            [0, 4.0, 5.0, 7.0],
            [4.0, 0, 7.0, 7.0],
            [5.0, 7.0, 0, 10.0],
            [7.0, 7.0, 10.0, 0]])
    form_objects = [
            Form.Matrix('matrix', 'distance matrix',
                D, MatrixUtil.assert_predistance)]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # read the distance matrix
    D = fs.matrix
    n = len(D)
    # get the list of implied variances
    V = [sum(row) / (n - 2) for row in D]
    # create the sigma matrix
    sigma = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            sigma[i][j] = (V[i] + V[j] - D[i][j]) / 2
    # begin the response
    out = StringIO()
    print >> out, MatrixUtil.m_to_string(sigma)
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
