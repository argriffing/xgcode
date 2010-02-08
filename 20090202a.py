"""Compare two matrices derived from the distance matrix.

This shows two ways of perturbing the inverse of a distance matrix
so that its rows and columns each end up summing to zero.
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
    # define the default distance matrix
    # this is from figure 2 of a paper called why neighbor joining works
    D = np.array([
        [0.0, 3.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0],
        [3.0, 0.0, 3.0, 3.0, 3.0, 2.0, 2.0, 2.0],
        [2.0, 3.0, 0.0, 0.1, 0.4, 3.0, 3.0, 3.0],
        [2.0, 3.0, 0.1, 0.0, 0.4, 3.0, 3.0, 3.0],
        [2.0, 3.0, 0.4, 0.4, 0.0, 3.0, 3.0, 3.0],
        [3.0, 2.0, 3.0, 3.0, 3.0, 0.0, 0.1, 0.4],
        [3.0, 2.0, 3.0, 3.0, 3.0, 0.1, 0.0, 0.4],
        [3.0, 2.0, 3.0, 3.0, 3.0, 0.4, 0.4, 0.0]])
    form_objects = [
            Form.Matrix('matrix', 'distance matrix',
                D, MatrixUtil.assert_predistance)]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # read the matrix
    D = fs.matrix
    n = len(D)
    if n < 3:
        raise HandlingError('the matrix should have at least three rows')
    # define the other matrices
    D_inv = np.linalg.inv(D)
    row_sums = np.sum(D_inv, 0)
    grand_sum = np.sum(D_inv)
    A = np.zeros((n,n))
    B = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            A[i][j] = row_sums[i] + row_sums[j] - grand_sum
            B[i][j] = row_sums[i] * row_sums[j] / grand_sum
    C = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            C[i][j] = D_inv[i][j] - B[i][j]
    # define the response
    out = StringIO()
    print >> out, 'additive:'
    print >> out, MatrixUtil.m_to_string(A)
    print >> out, 'multiplicative:'
    print >> out, MatrixUtil.m_to_string(B)
    for row in C:
        print >> out, sum(row)
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
