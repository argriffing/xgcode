"""Do a Procrustes transformation.

Get the Procrustes transformation of a list of input points
to match a list of reference points.
"""

from StringIO import StringIO
import math

import numpy as np

from SnippetUtil import HandlingError
import MatrixUtil
import Form

def todec(degrees, minutes):
    """
    @param degrees: latitude or longitude degrees
    @param minutes: latitude or longitude minutes
    @return: a floating point number
    """
    return degrees + minutes / 60.0

def get_form():
    """
    @return: the body of a form
    """
    # define some input points to transform
    input_points = np.array([
            (-4.5522481406, 10.4157358089),
            (-15.9814450114, -3.58894852685),
            (5.63165500818, 1.74036363828),
            (3.82188866923, -10.3582135754),
            (11.0801494746, 1.79106265508)])
    # define some reference points
    reference_points = np.array([
            (-todec(78, 39), todec(35, 46)),
            (-todec(84, 23), todec(33, 45)),
            (-todec(77, 02), todec(38, 53)),
            (-todec(79, 57), todec(40, 27)),
            (-todec(75, 10), todec(39, 57))])
    # define the form objects
    form_objects = [
            Form.Matrix('input_points', 'input points to transform',
                input_points),
            Form.Matrix('reference_points', 'reference points',
                reference_points)]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # do Procrustes Rotation as in section 14.7 of Multivariate Analysis by Mardia et al.
    Y_uncentered = fs.input_points
    X_uncentered = fs.reference_points
    if Y_uncentered.shape[0] != X_uncentered.shape[0]:
        raise HandlingError('the number of input points should be the same as the number of reference points')
    if Y_uncentered.shape[1] != X_uncentered.shape[1]:
        raise HandlingError('the dimensionality of the input points should be the same as the dimensionality of the reference points')
    n = len(Y_uncentered)
    # get the center of the reference points
    reference_center = sum(X_uncentered) / float(n)
    # center the input points and the reference points
    H = np.eye(n) - np.ones((n, n)) / n
    Y = np.dot(H, Y_uncentered)
    X = np.dot(H, X_uncentered)
    # get the singular value decomposition of a matrix
    Z = np.dot(Y.T, X)
    V, gamma, Uh = np.linalg.svd(Z)
    # get an orthogonal matrix that is part of the transformation
    A_hat = np.dot(V, Uh)
    # get a scaling factor
    c_hat = sum(gamma) / np.trace(np.dot(Y, Y.T))
    # begin the response
    out = StringIO()
    for point in Y:
        new_point = c_hat * np.dot(point, A_hat) + reference_center
        print >> out, '\t'.join(str(v) for v in new_point)
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
