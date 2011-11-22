"""
Check the eigendecomposition of a matrix constructed from a kernel function.

The kernel function is apparently the characteristic function
of the hyperbolic cosecant distribution.
It is f(t) = t / sinh(t).
The idea is that if you feed it a vector v then
you can get a matrix Mij = (vj - vi) / sinh(vj - vi)
such that M is guaranteed to be postive definite
as long as v was not something like zero.
"""

from StringIO import StringIO
import argparse
import math

import numpy as np
import scipy
from scipy import linalg

import Form
import FormOut
from MatrixUtil import ndot

#TODO write a rate matrix module now that I know how to use numpy

def get_form():
    """
    @return: the body of a form
    """
    # define the form objects
    form_objects = [
            Form.Integer('nstates', 'number of states', 4, low=2, high=12)]
    return form_objects

def get_form_out():
    return FormOut.Report()

def R_to_distn(R):
    """
    @param R: rate matrix
    @return: stationary distribution
    """
    n = len(R)
    Wl, Vl = scipy.linalg.eig(R, left=True, right=False)
    val_vec_pairs = [(abs(Wl[i]), Vl[:,i]) for i in range(n)]
    dummy, pi_eigenvector = min(val_vec_pairs)
    total = np.sum(pi_eigenvector)
    pi_arr = np.array([v/total for v in pi_eigenvector])
    return pi_arr

def sample_distribution(n):
    """
    @param n: number of states
    """
    # Get a nonnegative vector.
    v = np.random.rand(n)
    # Divide the vector by its nonnegative sum.
    distn = v / np.sum(v)
    return distn

def sample_symmetric_rate_matrix(n):
    """
    @param n: number of states
    """
    # Get a nonnegative asymmetric matrix.
    M = np.random.rand(n, n)
    # Symmetrize by adding to the transpose.
    S = M + M.T
    # Subtract row sum from diagonal.
    R = S - np.diag(np.sum(S, axis=1))
    return R

def to_gtr_c(S, v):
    """
    @param S: symmetric rate matrix
    @param v: target stationary distribution
    @return: a rate matrix with the target stationary distribution
    """
    # get the number of states
    n = len(v)
    # copy the symmetric rate matrix
    R = S.copy()
    # adjust the entries of the rate matrix
    for i in range(n):
        for j in range(n):
            R[i, j] *= math.sqrt(v[j] / v[i])
    # reset the diagonal entries of the rate matrix
    R -= np.diag(np.sum(R, axis=1))
    return R

def get_response_content(fs):
    np.set_printoptions(linewidth=200)
    out = StringIO()
    # do the analysis
    n = fs.nstates
    pi_m = sample_distribution(n)
    pi_q = sample_distribution(n)
    v = np.log(np.sqrt(pi_m / pi_q))
    K = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            x = v[j] - v[i]
            if x:
                K[i, j] = x / math.sinh(x)
            else:
                K[i, j] = 1.0
    W, V = scipy.linalg.eigh(K)
    # make a gtr mutation matrix
    S_precursor = sample_symmetric_rate_matrix(n)
    M = to_gtr_c(S_precursor, pi_m)
    M_distn = R_to_distn(M)
    if not np.allclose(M_distn, pi_m):
        raise ValueError('stationary distribution error')
    # resymmetrize
    lam = np.diag(np.sqrt(pi_m))
    rlam = np.diag(np.reciprocal(np.sqrt(pi_m)))
    S = ndot(lam, M, rlam)
    R = S * K
    lam = np.diag(np.sqrt(pi_q))
    rlam = np.diag(np.reciprocal(np.sqrt(pi_q)))
    Q_from_R = ndot(rlam, R, lam)
    Q_from_R -= np.diag(np.sum(Q_from_R, axis=1))
    Q_from_S = ndot(rlam, S, lam)
    Q_from_S -= np.diag(np.sum(Q_from_S, axis=1))
    Q_from_precursor = to_gtr_c(S_precursor, pi_q)
    # write the report
    print >> out, 'mutation process stationary distribution:'
    print >> out, pi_m
    print >> out
    print >> out, 'selection process stationary distribution:'
    print >> out, pi_q
    print >> out
    print >> out, 'vector to which the kernel function is applied:'
    print >> out, v
    print >> out
    print >> out, 'kernel matrix K:'
    print >> out, K
    print >> out
    print >> out, 'eigenvalues of K:'
    print >> out, W
    print >> out
    print >> out, 'eigenvectors of K:'
    print >> out, V
    print >> out
    print >> out, 'symmetric precursor matrix:'
    print >> out, S_precursor
    print >> out
    print >> out, 'rate matrix M:'
    print >> out, M
    print >> out
    print >> out, 'symmetrization S of rate matrix M:'
    print >> out, S
    print >> out
    print >> out
    print >> out, 'symmetrization R = S o K'
    print >> out, R
    print >> out
    print >> out, 'de-symmetrized rate matrix derived from R:'
    print >> out, Q_from_R
    print >> out
    print >> out
    print >> out, 'de-symmetrized rate matrix derived from S:'
    print >> out, Q_from_S
    print >> out
    print >> out
    print >> out, 'rate matrix derived from precursor rate matrix:'
    print >> out, Q_from_precursor
    print >> out
    return out.getvalue()

