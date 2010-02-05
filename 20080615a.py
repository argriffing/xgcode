"""Convert a path resistance matrix to an edge resistor matrix.
"""

from StringIO import StringIO

import numpy

from SnippetUtil import HandlingError
import MatrixUtil
import Euclid
import Form

def get_form():
    """
    @return: the body of a form
    """
    D = numpy.array([
        [0.0, 2.0, 2.0],
        [2.0, 0.0, 2.0],
        [2.0, 2.0, 0.0]])
    return [Form.Matrix('matrix', 'path resistance matrix', D, MatrixUtil.assert_predistance)]

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # read the distance matrix
    D = fs.matrix
    L = Euclid.edm_to_laplacian(D)
    resistor = -1/L
    resistor -= numpy.diag(numpy.diag(resistor))
    # write the edge resistor matrix
    out = StringIO()
    print >> out, MatrixUtil.m_to_string(resistor)
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
