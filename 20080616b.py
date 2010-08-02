"""Convert a path resistance matrix to a weighted adjacency matrix.
"""

import numpy as np

import MatrixUtil
import Euclid
import Form
import FormOut

def get_form():
    """
    @return: the body of a form
    """
    D = np.array([
            [0, 2, 2],
            [2, 0, 2],
            [2, 2, 0]])
    form_objects = [
            Form.Matrix('matrix', 'path resistance matrix',
                D, MatrixUtil.assert_predistance)]
    return form_objects

def get_form_out():
    return FormOut.Matrix()

def get_response_content(fs):
    D = fs.matrix
    A = Euclid.edm_to_adjacency(D)
    return MatrixUtil.m_to_string(A) + '\n'
