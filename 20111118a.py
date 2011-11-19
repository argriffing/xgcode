"""
Check ways to make a time-reversible rate matrix from a symmetric rate matrix.

The constructed rate matrix should have a specified stationary distribution,
and the symmetric rate matrix should be reconstructible from
the constructed rate matrix.
"""

from StringIO import StringIO
import argparse
import math

import numpy as np
import scipy
from scipy import linalg

import Form
import FormOut

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

def assert_reversible(R):
    if not is_reversible(R):
        msg = 'the rate matrix is not time-reversible'
        raise ValueError(msg)

def is_reversible(R):
    v = R_to_distn(R)
    n = len(v)
    for i in range(n):
        for j in range(i-1):
            forward = v[j] * R[j, i]
            backward = v[i] * R[i, j]
            if not np.allclose(forward, backward):
                return False
    return True

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

def sample_distribution(n):
    """
    @param n: number of states
    """
    # Get a nonnegative vector.
    v = np.random.rand(n)
    # Divide the vector by its nonnegative sum.
    distn = v / np.sum(v)
    return distn

def to_gtr_a(S, v):
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
            R[i, j] *= v[j]
    # reset the diagonal entries of the rate matrix
    R -= np.diag(np.sum(R, axis=1))
    return R

def to_gtr_b(S, v):
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
            R[i, j] /= v[i]
    # reset the diagonal entries of the rate matrix
    R -= np.diag(np.sum(R, axis=1))
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

def to_gtr_d(S, v):
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
    for a in range(n):
        for b in range(n):
            if a != b:
                # This equation is unnecessarily verbose
                # due to symmetry of S.
                # It should also work for asymmetric input rate matrices.
                tau = (v[b] * S[b, a]) / (v[a] * S[a, b])
                R[a, b] *= math.log(tau) / (1 - 1/tau)
    # reset the diagonal entries of the rate matrix
    R -= np.diag(np.sum(R, axis=1))
    return R

def to_gtr_e(S, v):
    """
    @param S: symmetric rate matrix
    @param v: target stationary distribution
    @return: a rate matrix with the target stationary distribution
    """
    p = R_to_distn(S)
    # get the number of states
    n = len(v)
    # copy the symmetric rate matrix
    R = S.copy()
    # adjust the entries of the rate matrix
    for a in range(n):
        for b in range(n):
            if a != b:
                # This equation is unnecessarily verbose
                # due to symmetry of S.
                # It should also work for asymmetric input rate matrices.
                tau = (v[b] / p[b]) / (v[a] / p[a])
                R[a, b] *= math.log(tau) / (1 - 1/tau)
    # reset the diagonal entries of the rate matrix
    R -= np.diag(np.sum(R, axis=1))
    return R


def to_s_a(R):
    """
    @param R: a time-reversible rate matrix
    """
    # get the number of states
    n = len(R)
    # get the stationary distribution of the rate matrix
    # TODO finish
    

def get_response_content(fs):
    np.set_printoptions(linewidth=200)
    out = StringIO()
    # do the analysis
    n = fs.nstates
    S = sample_symmetric_rate_matrix(n)
    v_uniform = np.ones(n, dtype=float) / n
    S_distn = R_to_distn(S)
    if not np.allclose(S_distn, v_uniform):
        msg = 'symmetric rate matrix should have uniform stationary distn'
        raise ValueError(msg)
    v = sample_distribution(n)
    # report the random inputs
    print >> out, 'random symmetric rate matrix:'
    print >> out, S
    print >> out
    print >> out, 'symmetric rate matrix stationary distribution:'
    print >> out, S_distn
    print >> out
    print >> out, 'random target stationary distribution:'
    print >> out, v
    print >> out
    print >> out
    for to_gtr in (to_gtr_a, to_gtr_b, to_gtr_c, to_gtr_d, to_gtr_e):
        R = to_gtr(S, v)
        R_distn = R_to_distn(R)
        if not np.allclose(R_distn, v):
            msg = 'constructed rate matrix should have target stationary distn'
            raise ValueError(msg)
        assert_reversible(R)
        # report the reconstruction
        print >> out, 'constructed rate matrix:'
        print >> out, R
        print >> out
        print >> out, 'constructed rate matrix stationary distribution:'
        print >> out, R_distn
        print >> out
        print >> out
    return out.getvalue()

