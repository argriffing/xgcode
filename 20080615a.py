"""Convert a path resistance matrix to an edge resistor matrix.
"""

import numpy as np

from SnippetUtil import HandlingError
import MatrixUtil
import Euclid
import Form
import FormOut

def get_form():
    """
    @return: the body of a form
    """
    D = np.array([
        [0.0, 2.0, 2.0],
        [2.0, 0.0, 2.0],
        [2.0, 2.0, 0.0]])
    form_objects = [
            Form.Matrix('matrix', 'path resistance matrix',
                D, MatrixUtil.assert_predistance)]
    return form_objects

def get_form_out():
    return FormOut.Matrix()

def get_response_content(fs):
    # read the distance matrix
    D = fs.matrix
    L = Euclid.edm_to_laplacian(D)
    resistor = -1/L
    resistor -= np.diag(np.diag(resistor))
    # return the edge resistor matrix
    return MatrixUtil.m_to_string(resistor) + '\n'
