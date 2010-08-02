"""Calculate the Moore-Penrose inverse of a matrix.
"""

import numpy as np

import Form
import FormOut
import MatrixUtil
import FelTree
import NewickIO
import const

g_default_string = const.read('20100730q')

def get_form():
    """
    @return: a list of form objects
    """
    tree = NewickIO.parse(g_default_string, FelTree.NewickTree)
    M = np.array(tree.get_full_distance_matrix())
    form_objects = [
            Form.Matrix('matrix', 'matrix', M),
            Form.Float('epsilon', 'show smaller elements as zero', '1e-10')]
    return form_objects

def get_form_out():
    return FormOut.Matrix()

def get_response_content(fs):
    # read the matrix
    M = fs.matrix
    # get the pseudo inverse
    M_pinv = np.linalg.pinv(M)
    # set small values to zero in the output
    M_pinv[abs(M_pinv) < fs.epsilon] = 0
    # return the response string
    return MatrixUtil.m_to_string(M_pinv) + '\n'
