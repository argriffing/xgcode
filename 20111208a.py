"""
Check spectra of specific rate matrices.
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
import mrate


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
    # define the exchangeability
    n = fs.nstates
    S = np.ones((n, n), dtype=float)
    S -= np.diag(np.sum(S, axis=1))
    # define the unnormalized target stationary distributions
    arr = []
    #
    v = np.ones(n)
    arr.append(v)
    #
    v = np.ones(n)
    v[0] += 0.1
    arr.append(v)
    #
    v = np.ones(n)
    v[0] -= 0.1
    arr.append(v)
    #
    v = np.ones(n)
    v[0] += 0.1
    v[1] += 0.1
    arr.append(v)
    #
    v = np.ones(n)
    v[0] -= 0.1
    v[1] -= 0.1
    #
    v = np.ones(n)
    v[0] += n
    arr.append(v)
    #
    v = np.ones(n)
    v[0] = 1.0 / n
    arr.append(v)
    #
    v = np.ones(n)
    v[0] += n
    v[1] += n
    arr.append(v)
    #
    v = np.ones(n)
    v[0] = 1.0 / n
    v[1] = 1.0 / n
    arr.append(v)
    #
    # write the report
    for v_weights in arr:
        v = v_weights / np.sum(v_weights)
        R = to_gtr_c(S, v)
        W, V = scipy.linalg.eig(R, left=True, right=False)
        recip_relaxation_time = sorted(abs(w) for w in W)[1]
        relaxation_time = 1.0 / recip_relaxation_time
        distn = mrate.R_to_distn(R)
        print >> out, 'observed stationary distribution:'
        print >> out, distn
        print >> out
        print >> out, 'relaxation time:'
        print >> out, relaxation_time
        print >> out
        print >> out, 'rate matrix:'
        print >> out, R
        print >> out
        print >> out, 'spectrum:'
        print >> out, W
        print >> out
        print >> out, 'left eigenvectors as columns:'
        print >> out, V
        print >> out
        print >> out
    return out.getvalue()

