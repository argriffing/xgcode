"""
Check a distinction between strain and stress in the context of MDS.

I think that the low dimensional Euclidean approximation found
by eigendecomposition of Gower's centered matrix minimizes strain
which has a less appealing interpretation than the minimization of stress.
Try some random things to see if there is really a difference
and if it means what I think it means.
"""

from StringIO import StringIO
import random

import numpy as np

import Form
import FormOut

def ndot(*args):
    M = args[0]
    for B in args[1:]:
        M = np.dot(M, B)
    return M

def double_centered(M):
    n = len(M)
    e = np.ones(n)
    I = np.eye(n)
    P = np.outer(e, e) / np.inner(e, e)
    H = I - P
    return ndot(H, M, H)

def get_random_points(npoints, ndim):
    """
    Each element will be a positive integer.
    Each row represents a point in Euclidean space.
    @param npoints: get this many points
    @param ndim: each point has this many dimensions
    @return: a float ndarray with shape (npoints, ndim)
    """
    X = np.zeros((npoints, ndim))
    for i in range(npoints):
        for j in range(ndim):
            X[i, j] = float(random.randrange(10))
    return X

def points_to_edm(X):
    """
    Get a Euclidean distance matrix.
    Entries are squared Euclidean distances between points.
    @param X: a float ndarray with shape (npoints, ndim)
    @return: a float ndarray with shape (npoints, npoints)
    """
    npoints, ndim = X.shape
    D = np.zeros((npoints, npoints))
    for i, vi in enumerate(X):
        for j, vj in enumerate(X):
            dv = vj - vi
            D[i, j] = np.dot(dv, dv)
    return D

def edm_to_gower_matrix(D):
    """
    @param D: an EDM
    @return: Gower's double centered matrix
    """
    G = -0.5*double_centered(D)
    return G


def get_form():
    """
    @return: the body of a form
    """
    # define the form objects
    form_objects = []
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):
    npoints = 5
    ndim = 4
    X = get_random_points(npoints, ndim)
    D = points_to_edm(X)
    G = edm_to_gower_matrix(D)
    # begin the output
    np.set_printoptions(linewidth=200)
    out = StringIO()
    print >> out, 'original points as rows:'
    print >> out, X
    print >> out, 'original EDM:'
    print >> out, D
    print >> out, 'original G:'
    print >> out, G
    return out.getvalue()

