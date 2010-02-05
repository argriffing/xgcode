"""Determine whether a Euclidean distance matrix is spherical.

A Euclidean distance matrix is spherical if and only if the sum of the elements of its pseudoinverse is positive.
This property is shown in
"Properties of Euclidean and Non-Euclidean distance matrices".
The example distance matrix is from the four taxon tree with unit branch lengths and with internal nodes.
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
        [0, 2, 3, 3, 1, 2],
        [2, 0, 3, 3, 1, 2],
        [3, 3, 0, 2, 2, 1],
        [3, 3, 2, 0, 2, 1],
        [1, 1, 2, 2, 0, 1],
        [2, 2, 1, 1, 1, 0]])
    # define the form objects
    return [Form.Matrix('matrix', 'Euclidean distance matrix', D, MatrixUtil.assert_predistance)]

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # read the matrix
    D = fs.matrix
    # begin the response
    out = StringIO()
    # look at the eigenvalues of the associated doubly centered covariance matrix
    HSH = Euclid.edm_to_dccov(D)
    w, V_T = numpy.linalg.eigh(HSH)
    V = V_T.T
    print >> out, 'eigenvalues of the associated doubly centered covariance matrix:'
    for x in reversed(sorted(w)):
        print >> out, x
    print >> out
    print >> out, 'eigenvector associated with last eigenvalue:'
    last_eigenvector = min(zip(w, V))[1]
    for x in last_eigenvector:
        print >> out, x
    print >> out
    # look at another criterion
    D_pinv = numpy.linalg.pinv(D)
    criterion = numpy.sum(D_pinv)
    if criterion > 0:
        print >> out, 'sum of elements of the pseudoinverse of the distance matrix is positive'
    else:
        print >> out, 'sum of elements of the pseudoinverse of the distance matrix is nonpositive'
    print >> out, 'A Euclidean distance matrix is spherical if and only if the sum of the elements of its pseudoinverse is positive.'
    print >> out, 'For this distance matrix, this sum is', criterion
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
