"""Do varimax rotation.
"""

import numpy as np

import MatrixUtil
import Form
import FormOut

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.Matrix('matrix', 'matrix',
                MatrixUtil.g_example_loading_matrix)]
    return form_objects

def get_form_out():
    return FormOut.Matrix()

def get_response_content(fs):
    rotated_matrix = np.dot(fs.matrix, MatrixUtil.varimax(fs.matrix))
    return MatrixUtil.m_to_string(rotated_matrix) + '\n'
