"""
Check spectra of related rate matrices.
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

def R_to_relaxation_time(R):
    n = len(R)
    Wl, Vl = scipy.linalg.eig(R, left=True, right=False)
    val_vec_pairs = [(abs(Wl[i]), Vl[:,i]) for i in range(n)]
    r_recip, fiedler = sorted(val_vec_pairs)[1]
    return 1 / r_recip

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

def symmetrized(R):
    """
    Get the symmetrized matrix.
    This returns a symmetric matrix that is not a rate matrix
    because rows do not sum to zero.
    """
    v = R_to_distn(R)
    lam = np.diag(np.sqrt(v))
    rlam = np.diag(np.reciprocal(np.sqrt(v)))
    return ndot(lam, R, rlam)

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
    # make a gtr mutation matrix
    L = sample_symmetric_rate_matrix(n)
    M = to_gtr_c(L, pi_m)
    M_sym = symmetrized(M)
    M_distn = R_to_distn(M)
    if not np.allclose(M_distn, pi_m):
        raise ValueError('mutation stationary distn error')
    Q = to_gtr_c(L, pi_q)
    Q_sym = symmetrized(Q)
    Q_distn = R_to_distn(Q)
    if not np.allclose(Q_distn, pi_q):
        raise ValueError('mutation-selection stationary distn error')
    # write the report
    print >> out, 'mutation process stationary distribution:'
    print >> out, pi_m
    print >> out
    print >> out, 'selection process stationary distribution:'
    print >> out, pi_q
    print >> out
    print >> out, 'symmetric rate matrix L:'
    print >> out, L
    print >> out
    print >> out, 'trace of L:'
    print >> out, np.trace(L)
    print >> out
    print >> out, 'eigendecomposition of L:'
    L_W, L_V = scipy.linalg.eig(L, right=False, left=True)
    print >> out, L_W
    print >> out, L_V
    print >> out
    print >> out, 'mutation rate matrix M:'
    print >> out, M
    print >> out
    print >> out, 'eigendecomposition of M:'
    M_W, M_V = scipy.linalg.eig(M, right=False, left=True)
    print >> out, M_W
    print >> out, M_V
    print >> out
    print >> out, '(symmetrized M) - (L):'
    print >> out, M_sym - L
    print >> out
    print >> out, 'trace of symmetrized M:'
    print >> out, np.trace(M_sym)
    print >> out
    print >> out, 'mutation-selection rate matrix Q:'
    print >> out, Q
    print >> out
    print >> out, 'eigendecomposition of Q:'
    Q_W, Q_V = scipy.linalg.eig(Q, right=False, left=True)
    print >> out, Q_W
    print >> out, Q_V
    print >> out
    print >> out, '(symmetrized Q) - (L):'
    print >> out, Q_sym - L
    print >> out
    print >> out, 'trace of symmetrized Q:'
    print >> out, np.trace(Q_sym)
    print >> out
    # look at how small changes in stationary distribution
    # affect the relaxation time
    eps = 1e-2
    G = L
    G_sym = L
    G_distn = np.ones(n) / float(n)
    r_mutation = R_to_relaxation_time(G)
    n = len(G)
    W, V = scipy.linalg.eigh(G_sym)
    val_vec_pairs = [(abs(W[i]), V[:,i]) for i in range(n)]
    for val, vec in sorted(val_vec_pairs):
        # get the target stationary distribution
        v = G_distn + vec*eps
        v /= sum(v)
        # get the mutation-selection rate matrix
        Q = to_gtr_c(L, v)
        # get the relaxation time
        r = R_to_relaxation_time(Q)
        stationary_delta = v - G_distn
        stationary_delta_mag = np.linalg.norm(stationary_delta)
        relaxation_time_derivative = (r - r_mutation) / stationary_delta_mag
        print >> out, 'eigenvalue:'
        print >> out, val
        print >> out, 'eigenvector:'
        print >> out, vec
        print >> out, 'stationary distribution delta:'
        print >> out, stationary_delta
        print >> out, 'relaxation time change per stationary distn change:'
        print >> out, relaxation_time_derivative
        print >> out
    print >> out
    print >> out
    for val, vec in sorted(val_vec_pairs):
        # get the target stationary distribution
        dvec = vec*vec - np.mean(vec*vec)
        #v = M_distn + dvec*eps
        v = np.exp(np.log(G_distn) + dvec*eps)
        v /= sum(v)
        # get the mutation-selection rate matrix
        Q = to_gtr_c(L, v)
        # get the relaxation time
        r = R_to_relaxation_time(Q)
        stationary_delta = v - G_distn
        stationary_delta_mag = np.linalg.norm(stationary_delta)
        #relaxation_time_derivative = (r - r_mutation) / stationary_delta_mag
        relaxation_time_log = math.log(r / r_mutation) / stationary_delta_mag
        print >> out, 'eigenvalue:'
        print >> out, val
        #print >> out, 'M eigenvector:'
        #print >> out, vec
        print >> out, 'centered entrywise square of eigenvector:'
        print >> out, dvec
        print >> out, 'stationary distribution delta:'
        print >> out, stationary_delta
        #print >> out, 'relaxation time change:'
        #print >> out, r, '-', r_mutation, '=', r - r_mutation
        #print >> out, 'relaxation time change per stationary distn change:'
        #print >> out, relaxation_time_derivative
        print >> out, 'log relaxation time ratio per stationary distn change:'
        print >> out, relaxation_time_log
        print >> out
    return out.getvalue()

