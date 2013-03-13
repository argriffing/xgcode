"""
Check a connection between Laplacian matrix and neighbor joining criteria.
"""

from StringIO import StringIO

import numpy as np
import scipy

import Form
import FormOut
import MatrixUtil
from MatrixUtil import ndot

def get_form():
    """
    @return: the body of a form
    """
    return [
            ]

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):

    # set up print options
    np.set_printoptions(
            linewidth=1000000,
            threshold=1000000,
            )
    out = StringIO()

    # define the Laplacian matrix
    L = np.array([
        [ 1,  0,  0,  0,  0, -1,  0,  0],
        [ 0,  2,  0,  0,  0, -2,  0,  0],
        [ 0,  0,  3,  0,  0,  0, -3,  0],
        [ 0,  0,  0,  2,  0,  0, -2,  0],
        [ 0,  0,  0,  0,  1,  0,  0, -1],
        [-1, -2,  0,  0,  0,  4,  0, -1],
        [ 0,  0, -3, -2,  0,  0,  6, -1],
        [ 0,  0,  0,  0, -1, -1, -1,  3],
        ], dtype=float)

    # remove the last two columns by schur complementation
    L_schur = L[:-2, :-2] - ndot(
            L[:-2, -2:], 
            scipy.linalg.inv(L[-2:, -2:]),
            L[-2:, :-2])

    # get the trailing block of the matrix
    L_schur_component = L_schur[-4:, -4:]

    # get the part corresponding to the inverse of a rooted covariance matrix
    L_schur_rooted = L_schur_component[:-1, :-1]
    
    # get the corresponding covariance matrix
    cov = scipy.linalg.inv(L_schur_rooted)

    # print the matrices
    print >> out, 'L:'
    print >> out, L
    print >> out
    print >> out, 'schur complement of two internal vertices (7, 8) in L:'
    print >> out, L_schur
    print >> out
    print >> out, 'a component (3, 4, 5, 6) of the schur complement:'
    print >> out, L_schur_component
    print >> out
    print >> out, 'a piece (3, 4, 5) of the component:'
    print >> out, L_schur_rooted
    print >> out
    print >> out, 'the corresponding rooted covariance matrix:'
    print >> out, cov
    print >> out
    print >> out, 'trace of covariance matrix:'
    print >> out, np.trace(cov)
    print >> out

    # show the result
    return out.getvalue()

