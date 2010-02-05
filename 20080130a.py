"""Given a rate matrix and a time, compute the transition matrix.

The matrix exponential computation uses the default scipy implementation,
which I think is Pade approximation of order 7.
"""

from StringIO import StringIO

# the scipy version of linalg is necessary for expm
from scipy import linalg
import numpy

from SnippetUtil import HandlingError
import Newick
import MatrixUtil
import RateMatrix
import Form

def get_form():
    """
    @return: the body of a form
    """
    # define the default rate matrix
    R = numpy.array([
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

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # exponentiate the matrix
    R = fs.matrix * fs.time
    T = linalg.expm(R)
    # write the response
    out = StringIO()
    print >> out, MatrixUtil.m_to_string(T)
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
