"""
Check a distinction between strain and stress in the context of MDS.

I think that the low dimensional Euclidean approximation found
by eigendecomposition of Gower's centered matrix minimizes strain
which has a less appealing interpretation than the minimization of stress.
Try some random things to see if there is really a difference
and if it means what I think it means.
"strain" is the sum of squares of errors of elements
of Gower's centered matrix.
"stress" is the sum of squares of errors of square roots of elements
of the Euclidean distance matrix.
"stress b" is the sum of squares of errors of elements
of the Euclidean distance matrix.
Note that elements of a "Euclidean distance matrix" (EDM) are
squares of Euclidean distances.
"""

from StringIO import StringIO
import random

import numpy as np
import scipy
from scipy import linalg
from scipy import optimize

import Form
import FormOut
from MatrixUtil import ndot
from MatrixUtil import double_centered_slow

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
    G = -0.5 * double_centered_slow(D)
    return G

def gower_matrices_to_strain(Ga, Gb):
    """
    Strain is the sum of squares of Gower matrix errors.
    @param Ga: a Gower matrix
    @param Gb: another Gower matrix
    @return: the sum of squares of Gower matrix errors
    """
    dG = Ga - Gb
    strain = np.sum(dG * dG)
    return strain

def edm_matrices_to_stress(Da, Db):
    """
    Stress is the sum of squares of Euclidean distance errors.
    From Wikipedia,
    stress(X) = sum_{i<j<=n}{(d_ij(X) - del_ij)^2
    where d_ij(X) is the Euclidean distance between i and j.
    @param Da: an EDM whose entries are squares of Euclidean distances
    @param Db: another EDM whose entries are squares of Euclidean distances
    """
    M = np.sqrt(Db) - np.sqrt(Da)
    stress = np.sum(M * M)
    return stress

def edm_matrices_to_stress_b(Da, Db):
    """
    Use a different definition of stress.
    This is the sum of squares of errors of squared Euclidean distances.
    @param Da: an EDM whose entries are squares of Euclidean distances
    @param Db: another EDM whose entries are squares of Euclidean distances
    """
    M = Db - Da
    stress_b = np.sum(M * M)
    return stress_b

def get_form():
    """
    @return: the body of a form
    """
    # define the form objects
    form_objects = []
    return form_objects

def get_form_out():
    return FormOut.Report()

class CurriedStress:
    def __init__(self, D, X_shape):
        """
        @param D: EDM
        @param X_shape: the shape of the point matrix
        """
        self.D = D
        self.X_shape = X_shape
    def __call__(self, X_flat):
        """
        @param X_flat: flattened 2D array
        @return: error
        """
        X = np.reshape(X_flat, self.X_shape)
        D = points_to_edm(X)
        error = edm_matrices_to_stress(self.D, D)
        return error

class CurriedStressB:
    def __init__(self, D, X_shape):
        """
        @param D: EDM
        @param X_shape: the shape of the point matrix
        """
        self.D = D
        self.X_shape = X_shape
    def __call__(self, X_flat):
        """
        @param X_flat: flattened 2D array
        @return: error
        """
        X = np.reshape(X_flat, self.X_shape)
        D = points_to_edm(X)
        error = edm_matrices_to_stress_b(self.D, D)
        return error

class CurriedStrain:
    def __init__(self, G, X_shape):
        """
        @param G: Gower matrix
        @param X_shape: the shape of the point matrix
        """
        self.G = G
        self.X_shape = X_shape
    def __call__(self, X_flat):
        """
        @param X_flat: flattened 2D array
        @return: error
        """
        X = np.reshape(X_flat, self.X_shape)
        D = points_to_edm(X)
        G = edm_to_gower_matrix(D)
        error = gower_matrices_to_strain(self.G, G)
        return error

def get_summary_string(X_true, X_test):
    # define the true values
    D_true = points_to_edm(X_true)
    G_true = edm_to_gower_matrix(D_true)
    # define the test values
    D = points_to_edm(X_test)
    G = edm_to_gower_matrix(D)
    # get the errors
    strain = gower_matrices_to_strain(G, G_true)
    stress = edm_matrices_to_stress(D, D_true)
    stress_b = edm_matrices_to_stress_b(D, D_true)
    # get the summary string
    out = StringIO()
    print >> out, 'points as rows:'
    print >> out, X_test
    print >> out, 'EDM:'
    print >> out, D
    print >> out, 'G:'
    print >> out, G
    print >> out, 'strain:'
    print >> out, strain
    print >> out, 'stress:'
    print >> out, stress
    print >> out, 'stress b:'
    print >> out, stress_b
    return out.getvalue().rstrip()

def get_response_content(fs):
    np.set_printoptions(linewidth=200)
    npoints = 5
    ndim = 4
    r = 2
    # get the original points
    X = get_random_points(npoints, ndim)
    D = points_to_edm(X)
    G = edm_to_gower_matrix(D)
    # get a low dimensional approximation
    W, V = scipy.linalg.eigh(G)
    V_r = V[:, -r:]
    W_r = W[-r:]
    X_r = np.dot(V_r, np.diag(np.sqrt(W_r)))
    # numerically minimize strain
    X_flat = np.random.rand(npoints*r)
    f = CurriedStrain(G, (npoints, r))
    result = scipy.optimize.fmin(f, X_flat, disp=0, full_output=1)
    value, error, niter, ncalls, warn = result
    X_min_strain = value.reshape((npoints, r))
    # numerically minimize stress
    X_flat = np.random.rand(npoints*r)
    f = CurriedStress(D, (npoints, r))
    result = scipy.optimize.fmin(f, X_flat, disp=0, full_output=1)
    value, error, niter, ncalls, warn = result
    X_min_stress = value.reshape((npoints, r))
    # numerically minimize stress b
    X_flat = np.random.rand(npoints*r)
    f = CurriedStressB(D, (npoints, r))
    result = scipy.optimize.fmin(f, X_flat, disp=0, full_output=1)
    value, error, niter, ncalls, warn = result
    X_min_stress_b = value.reshape((npoints, r))
    # begin the output
    out = StringIO()
    print >> out, 'original:'
    print >> out, get_summary_string(X, X)
    print >> out
    print >> out, 'closed form low dimensional approximation:'
    print >> out, get_summary_string(X, X_r)
    print >> out
    print >> out, 'stochastic min strain low dimensional approximation:'
    print >> out, get_summary_string(X, X_min_strain)
    print >> out
    print >> out, 'stochastic min stress low dimensional approximation:'
    print >> out, get_summary_string(X, X_min_stress)
    print >> out
    print >> out, 'stochastic min stress b low dimensional approximation:'
    print >> out, get_summary_string(X, X_min_stress_b)
    return out.getvalue()

