"""Given a rate matrix and a time, compute the transition matrix.

The matrix exponential computation uses the default scipy implementation,
which I think is Pade approximation of order 7.
"""

from scipy import linalg
import numpy as np

from SnippetUtil import HandlingError
import Newick
import MatrixUtil
import RateMatrix
import Form
import FormOut

def get_form():
    """
    @return: the body of a form
    """
    # define the default rate matrix
    R = np.array([
        [-1, 1/3.0, 1/3.0, 1/3.0],
        [1/3.0, -1, 1/3.0, 1/3.0],
        [1/3.0, 1/3.0, -1, 1/3.0],
        [1/3.0, 1/3.0, 1/3.0, -1]])
    # define the form objects
    form_objects = [
            Form.Matrix('matrix', 'rate matrix', 
                R, MatrixUtil.assert_rate_matrix),
            Form.Float('time', 'time', 2.0, low_exclusive=0)]
    return form_objects

def get_form_out():
    return FormOut.TransitionMatrix()

def get_response_content(fs):
    R = fs.matrix * fs.time
    T = linalg.expm(R)
    return MatrixUtil.m_to_string(T) + '\n'
