"""Seek a matrix where index merging and Schur complementation do not commute.

Look for a matrix where index merging and Schur complementation do not commute.
"""

# Matrices in this script are implemented as two dimensional numpy arrays.

from StringIO import StringIO
import random

import numpy as np

from SnippetUtil import HandlingError
import Form
import SchurAlgebra

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.Integer('block_size', 'block size',
                2, low=1, high=6)]
    return form_objects

def sample_matrix(block_size):
    """
    Sample a block 4x4 symmetric matrix.
    Each block is square with block_size rows.
    @param block_size: the number of rows in square blocks of the partitioned matrix
    @return: a sampled matrix
    """
    n = block_size * 4
    expected_value = 1.0
    M = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            M[i,j] = random.expovariate(expected_value)
    return M

def analyze_matrix(M, block_size):
    """
    Get results for g(f(M)) and f(g(M)).
    @param M: a matrix
    @param block_size: the number of rows in square blocks of the partitioned matrix
    @return: a string of results
    """
    # define the response
    out = StringIO()
    # get the new matrix using the first composition of functions
    M_11 = SchurAlgebra.mmerge(M, set(range(2*block_size)))
    M_12 = SchurAlgebra.mschur(M_11, set(1 + block_size + k for k in range(block_size)))
    print >> out, M_12
    # get the new matrix using the second composition of functions
    M_21 = SchurAlgebra.mschur(M, set(3*block_size + k for k in range(block_size)))
    M_22 = SchurAlgebra.mmerge(M_21, set(range(2*block_size)))
    print >> out, M_22
    if np.allclose(M_12, M_22):
        print >> out, 'the matrices are similar'
    else:
        print >> out, 'the matrices are different'
    return out.getvalue().strip()

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # define the response
    out = StringIO()
    # sample a matrix
    b = fs.block_size
    M = sample_matrix(b)
    print >> out, 'asymmetric analysis:'
    print >> out, analyze_matrix(M, b)
    print >> out
    print >> out, 'symmetric analysis:'
    print >> out, analyze_matrix(M + M.T, b)
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
