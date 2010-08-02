"""Normalize a rate matrix to have one expected transition per unit time.
"""

import numpy as np

from SnippetUtil import HandlingError
import MatrixUtil
import RateMatrix
import Form
import FormOut

#FIXME numpy may not be necessary

def get_form():
    """
    @return: the body of a form
    """
    # define the default matrix string
    R = np.array([
            [-3, 1, 1, 1],
            [1, -3, 1, 1],
            [1, 1, -3, 1],
            [1, 1, 1, -3]])
    # define the form objects
    form_objects = [
            Form.Matrix('matrix', 'rate matrix',
                R, MatrixUtil.assert_rate_matrix)]
    return form_objects

def get_form_out():
    return FormOut.RateMatrix()

def get_response_content(fs):
    # read the matrix from the form data
    R = fs.matrix
    n = len(R)
    # convert the row major rate matrix to a rate matrix object
    arbitrary_states = [str(x) for x in range(n)]
    rate_matrix_object = RateMatrix.RateMatrix(R.tolist(), arbitrary_states)
    rate_matrix_object.normalize()
    normalized_row_major = rate_matrix_object.get_row_major_rate_matrix()
    # return the rate matrix
    return MatrixUtil.m_to_string(normalized_row_major) + '\n'
