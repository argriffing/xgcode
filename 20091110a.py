"""Do stress majorization. [UNFINISHED]

Start by finding the 2D MDS of a Laplacian matrix.
Then do a few iterations of the SMACOF algorithm
in an attempt to reduce the stress objective function.
Finally, compute for the MDS solution and for the stress-minimized solution.
The stress-minimized solution should have less stress.
"""


from StringIO import StringIO
import math

import numpy as np

from SnippetUtil import HandlingError
import SnippetUtil
import Form
import FormOut
import NewickIO
import NeighborJoining
import FelTree
import Steiner

g_laplacian = np.array([
    [ 1.0, -1.0,  0.0,  0.0],
    [-1.0,  2.0, -1.0,  0.0],
    [ 0.0, -1.0,  2.0, -1.0],
    [ 0.0,  0.0, -1.0,  1.0]])

def get_form():
    """
    @return: a list of form objects
    """
    # define the list of form objects
    form_objects = [
            Form.Matrix('laplacian', 'laplacian matrix', g_laplacian)]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_stress(D_expected, D_observed):
    """
    Stress is a monotonic function of the Frobenius norm of the difference of EDMs.
    """
    if D_expected.shape != D_observed.shape:
        raise ValueError('the distance matrices should have the same shape')
    n = len(D_expected)
    stress = 0
    for i in range(n):
        for j in range(i):
            stress += (D_expected[i,j] - D_observed[i,j])**2
    return stress

def do_smacof(D):
    """
    @param D: a matrix of squared Euclidean distances
    @return: a matrix of 2D points
    """
    pass

def get_response_content(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    return 'not implemented\n'
