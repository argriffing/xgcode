"""
Check the eigendecomposition of a parent independent Markov process.
"""

from StringIO import StringIO
import random
import math

import numpy as np
import scipy
from scipy import linalg

import Form
import FormOut
import mrate
import ctmcmi
import combobreaker
import MatrixUtil
from MatrixUtil import ndot


def get_form():
    form_objects = [
            Form.Integer('nstates', 'number of states', 4, low=2, high=10)]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_mi(Q, v, t):
    """
    Try to get a closed form for the mutual information.
    @param Q: parent independent rate matrix
    @param v: stationary distribution
    @param t: amount of time
    """
    n = len(v)
    mi = 0
    for i in range(n):
        # get probability of randomization
        rate = -Q[i,i]
        randomization_rate = rate / (1 - v[i])
        pnr = math.exp(-randomization_rate*t)
        pr = 1 - pnr
        # add the mutual information contribution of no net change
        pii = pnr + (1 - pnr) * v[i]
        mi += v[i] * pii * math.log(pii / v[i])
        # add the mutual information contribution of net change
        mi += v[i] * (1 - v[i]) * pr * math.log(pr)
    return mi

def get_transition_matrix(Q, v, t):
    """
    Try to get a closed form for the transition probability matrix.
    @param Q: parent independent rate matrix
    @param v: stationary distribution
    @param t: amount of time
    """
    n = len(v)
    P = np.zeros_like(Q)
    for i in range(n):
        rate = -Q[i,i]
        randomization_rate = rate / (1 - v[i])
        # compute the probability of no randomization before time t
        p_no_random = math.exp(-randomization_rate*t)
        # if there is no randomization we stay where we are
        P[i, i] += p_no_random
        # if there is randomization we could end up anywhere
        P[i] += (1 - p_no_random) * v
    return P

def process(nstates):
    np.set_printoptions(linewidth=200)
    # Sample a rate matrix.
    # Use a trick by Robert Kern to left and right multiply by diagonals.
    # http://mail.scipy.org/pipermail/numpy-discussion/2007-March/
    # 026809.html
    S = MatrixUtil.sample_pos_sym_matrix(nstates)
    v = mrate.sample_distn(nstates)
    R = (v**-0.5)[:,np.newaxis] * S * (v**0.5)
    R -= np.diag(np.sum(R, axis=1))
    # Construct a parent-independent process
    # with the same max rate and stationary distribution
    # as the sampled process.
    Q = np.outer(np.ones(nstates), v)
    Q -= np.diag(np.sum(Q, axis=1))
    rescaling_factor = max(np.diag(R) / np.diag(Q))
    Q *= rescaling_factor
    # sample a random time
    time_mu = 0.01
    t = random.expovariate(1 / time_mu)
    # Check that the mutual information of the
    # parent independent process is smaller.
    out = StringIO()
    print >> out, 'sampled symmetric part of the rate matrix S:'
    print >> out, S
    print >> out
    print >> out, 'sampled stationary distribution v:'
    print >> out, v
    print >> out
    print >> out, 'sqrt stationary distribution:'
    print >> out, np.sqrt(v)
    print >> out
    print >> out, 'implied rate matrix R:'
    print >> out, R
    print >> out
    print >> out, 'related parent-independent rate matrix Q:'
    print >> out, Q
    print >> out
    print >> out, 'rescaling factor:'
    print >> out, rescaling_factor
    print >> out
    print >> out, 'eigenvalues of Q:'
    print >> out, scipy.linalg.eigvals(Q)
    print >> out
    print >> out, 'symmetric matrix similar to Q:'
    S = ndot(np.diag(np.sqrt(v)), Q, np.diag(1/np.sqrt(v)))
    print >> out, S
    print >> out
    print >> out, 'eigendecomposition of the similar matrix:'
    W, V = scipy.linalg.eigh(S)
    print >> out, V
    print >> out, np.diag(W)
    print >> out, V.T
    print >> out
    #
    print >> out, 'time:', t
    print >> out
    print >> out, 'stationary distribution entropy:', -ndot(v, np.log(v))
    print >> out
    # 
    P_by_hand = get_transition_matrix(Q, v, t)
    print >> out, 'parent independent transition matrix computed by hand:'
    print >> out, P_by_hand
    print >> out
    print >> out, 'parent independent transition matrix computed by expm:'
    print >> out, scipy.linalg.expm(Q*t)
    print >> out
    #
    print >> out, 'parent independent m.i. by hand:'
    print >> out, get_mi(Q, v, t)
    print >> out
    print >> out, 'parent independent m.i. by expm:'
    print >> out, ctmcmi.get_expected_ll_ratio(Q, t)
    print >> out
    #
    print >> out, 'original process m.i. by expm:'
    print >> out, ctmcmi.get_expected_ll_ratio(R, t)
    print >> out
    #
    return out.getvalue().rstrip()

def get_response_content(fs):
    return process(fs.nstates) + '\n'

