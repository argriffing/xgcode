"""Convert a path resistance matrix to a weighted adjacency matrix.
"""

from StringIO import StringIO

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
            [0, 2, 2],
            [2, 0, 2],
            [2, 2, 0]])
    form_objects = [
            Form.Matrix('matrix', 'path resistance matrix',
                D, MatrixUtil.assert_predistance)]
    return form_objects

def get_form_out():
    return FormOut.Matrix()

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # read the matrix
    D = fs.matrix
    # create the weighted adjacency matrix
    A = Euclid.edm_to_adjacency(D)
    # start to prepare the reponse
    out = StringIO()
    print >> out, MatrixUtil.m_to_string(A)
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
