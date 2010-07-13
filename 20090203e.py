"""Given a list of points, find the Euclidean distance matrix.

Note that the elements of the Euclidean distance matrix (EDM) are
usually defined as squared pairwise Euclidean distances.
"""

from StringIO import StringIO
import math

import numpy as np

from SnippetUtil import HandlingError
import MatrixUtil
from Form import RadioItem
import Form
import FormOut

def get_form():
    """
    @return: the body of a form
    """
    M = np.array([
            [0.701808823709, -0.66990306947, 0.226340180955, 0.166666666667, -0.0863966143067],
            [0.701808823709, 0.66990306947, -0.226340180955, 0.166666666667, -0.0863966143067],
            [-0.701808823709, 0.226340180955, 0.66990306947, 0.166666666667, 0.0863966143067],
            [-0.701808823709, -0.226340180955, -0.66990306947, 0.166666666667, 0.0863966143067],
            [0.394102719008, 0, 0, -0.333333333333, 0.307706104701],
            [-0.394102719008, 0, 0, -0.333333333333, -0.307706104701]])
    # define the form objects
    form_objects = [
            Form.Matrix('points', 'one point on each row', M),
            Form.RadioGroup('options', 'distance options', [
                Form.RadioItem('squared',
                    'get squared euclidean distances', True),
                Form.RadioItem('not_squared',
                    'get euclidean distances')])]
    return form_objects

def get_form_out():
    return FormOut.Matrix()

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # read the points
    M = fs.points
    # get the squared distance between each pair of points
    n = len(M)
    D = np.zeros((n,n))
    for i, pa in enumerate(M):
        for j, pb in enumerate(M):
            d = np.linalg.norm(pb - pa)
            if fs.squared:
                D[i][j] = d*d
            else:
                D[i][j] = d
    # define the response
    out = StringIO()
    print >> out, MatrixUtil.m_to_string(D)
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
