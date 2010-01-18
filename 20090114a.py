"""Check the three-point and four-point conditions for a distance matrix.
"""

import StringIO

import numpy

from SnippetUtil import HandlingError
import MatrixUtil
import Form

def get_form():
    """
    @return: the body of a form
    """
    D = numpy.array([
        [0, 4, 5, 7],
        [4, 0, 7, 7],
        [5, 7, 0, 10],
        [7, 7, 10, 0]])
    return [Form.Matrix('matrix', 'distance matrix', D, MatrixUtil.assert_predistance)]

def check_three_point_condition(D):
    """
    This is also known as the triangle inequality.
    @param D: a distance matrix
    @return: True if the three point condition is satisfied
    """
    n = len(D)
    for i in range(n):
        for j in range(n):
            for k in range(n):
                condition = D[i][k] <= D[i][j] + D[j][k]
                if not condition:
                    return False
    return True

def check_four_point_condition(D):
    """
    This is also known as metric additivity.
    @param D: a distance matrix
    @return: True if the three point condition is satisfied
    """
    n = len(D)
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for l in range(n):
                    condition = D[i][j] + D[k][l] <= max(D[i][k] + D[j][l], D[j][k] + D[i][l])
                    if not condition:
                        return False
    return True

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # read the matrix
    D = fs.matrix
    # begin the response
    out = StringIO.StringIO()
    print >> out, 'satisifies the three-point condition:', check_three_point_condition(D)
    print >> out, 'satisifies the four-point condition:', check_four_point_condition(D)
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
