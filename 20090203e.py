"""Given a list of points, find the Euclidean distance matrix.

Note that the elements of the Euclidean distance matrix (EDM) are
usually defined as squared pairwise Euclidean distances.
"""

import StringIO
import math

import numpy

from SnippetUtil import HandlingError
import MatrixUtil
import Form

def get_form():
    """
    @return: the body of a form
    """
    M = numpy.array([
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
                Form.RadioItem('squared', 'get the squared euclidean distances', True),
                Form.RadioItem('not_squared', 'get the euclidean distances')])]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # read the points
    M = fs.points
    # get the squared distance between each pair of points
    n = len(M)
    D = numpy.zeros((n,n))
    for i, pa in enumerate(M):
        for j, pb in enumerate(M):
            d = numpy.linalg.norm(pb - pa)
            if fs.squared:
                D[i][j] = d*d
            else:
                D[i][j] = d
    # define the response
    out = StringIO.StringIO()
    print >> out, MatrixUtil.m_to_string(D)
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
