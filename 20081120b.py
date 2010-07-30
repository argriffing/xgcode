"""Calculate the Moore-Penrose inverse of a matrix.
"""

from StringIO import StringIO

import numpy as np

from SnippetUtil import HandlingError
import SnippetUtil
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

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # read the matrix
    M = fs.matrix
    # get the pseudo inverse
    M_pinv = np.linalg.pinv(M)
    # set small values to zero in the output
    M_pinv[abs(M_pinv) < fs.epsilon] = 0
    # create the response string
    out = StringIO()
    print >> out, MatrixUtil.m_to_string(M_pinv)
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
