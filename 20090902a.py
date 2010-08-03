"""Look at properties of (D^-1 - L)^-1 in light of the NJ Q matrix.

The neighbor joining criterion involves a matrix that is like
L^+ + D/n, which is like a distance matrix perturbation of
the inverse of the Laplacian.
On the other hand Bapat et al. look at the inverse of
a matrix that is an inverse distance matrix perturbation of the laplacian.
Maybe there is some connection between these matrices, or maybe not.
"""

import numpy as np

from SnippetUtil import HandlingError
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
    # return the response
    return MatrixUtil.m_to_string(sigma) + '\n'
