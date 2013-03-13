"""
Given a distance matrix, get a covariance-like neighbor joining matrix.
"""

from StringIO import StringIO

import numpy as np
import scipy

import MatrixUtil
import Form
import FormOut

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

def get_form_out():
    return FormOut.Matrix()

def get_response_content(fs):
    out = StringIO()
    # read the distance matrix
    D = fs.matrix
    n = len(D)
    # get the list of implied variances
    V = [sum(row) / (n - 2) for row in D]
    # create the sigma matrix
    sigma = np.empty((n,n))
    for i in range(n):
        for j in range(n):
            sigma[i][j] = (V[i] + V[j] - D[i][j]) / 2.0
    # compute the eigendecomposition
    w, v = scipy.linalg.eigh(sigma)
    # write the response
    print >> out, 'covariance-like matrix:'
    print >> out, MatrixUtil.m_to_string(sigma)
    print >> out
    print >> out, 'eigenvalues:'
    print >> out, w
    print >> out
    print >> out, 'eigenvectors:'
    print >> out, v
    print >> out
    # return the response
    return out.getvalue().rstrip()
