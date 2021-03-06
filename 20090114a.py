"""Check the three-point and four-point conditions for a distance matrix.
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
    D = np.array([
        [0, 4, 5, 7],
        [4, 0, 7, 7],
        [5, 7, 0, 10],
        [7, 7, 10, 0]])
    form_objects = [
            Form.Matrix('matrix', 'distance matrix',
                D, MatrixUtil.assert_predistance)]
    return form_objects

def get_form_out():
    return FormOut.Report()

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

def get_response_content(fs):
    # read the matrix
    D = fs.matrix
    # begin the response
    out = StringIO()
    three_point = check_three_point_condition(D)
    four_point = check_four_point_condition(D)
    print >> out, 'satisifies the three-point condition:', three_point
    print >> out, 'satisifies the four-point condition:', four_point
    # write the response
    return out.getvalue()
