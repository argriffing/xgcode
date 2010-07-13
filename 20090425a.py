"""Do varimax rotation.
"""

from StringIO import StringIO

import numpy as np

from SnippetUtil import HandlingError
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

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    rotated_matrix = np.dot(fs.matrix, MatrixUtil.varimax(fs.matrix))
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, MatrixUtil.m_to_string(rotated_matrix)
