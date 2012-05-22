"""
Check the properties of simplifications of GTR Markov processes.

The two simplifications each preserve the stationary distribution
and a rate summary.
The 'parent-independent' simplification preserves the 'randomization rate.'
The 'child-independent' simplification preserves the expected rate.
"""

from StringIO import StringIO
import random
import math
import itertools
from itertools import combinations
from itertools import product

import numpy as np
import scipy
from scipy import linalg

import Form
import FormOut
import mrate
import ctmcmi
import cheeger
import msimpl
import iterutils
import MatrixUtil
from MatrixUtil import ndot


def get_form():
    form_objects = [
            Form.RadioGroup('simplification', 'simplification', [
                Form.RadioItem('parent_indep', 'parent independent', True),
                Form.RadioItem('bipartitioned', 'bipartitioned'),
                Form.RadioItem('child_indep', 'child independent')]),
            Form.Integer('nstates', 'number of states', 4, low=2, high=10),
            Form.Float('t', 'time', '0.1', low_exclusive=0)]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_randomization_rate(Q, v):
    """
    @param Q: reversible rate matrix
    @param v: stationary distribution
    @return: randomization rate
    """
    nstates = len(v)
    return max(get_randomization_candidate(Q, v, i) for i in range(nstates))

def get_randomization_candidate(Q, v, i):
    """
    @param Q: reversible rate matrix
    @param v: stationary distribution
    @param i: a state
    @return: the randomization rate for the current state
    """
    return -Q[i,i] / (1 - v[i])

def get_cheeger_constant(R, v):
    """
    This is also known as the second isoperimetric constant.
    @param R: a reversible rate matrix
    @param v: stationary distribution
    @return: the second isoperimetric constant
    """
    n = len(v)
    I2 = None
    for A_tuple in iterutils.powerset(range(n)):
        # define the vertex set and its complement
        A = set(A_tuple)
        B = set(range(n)) - A
        A_measure = sum(v[i] for i in A)
        B_measure = sum(v[i] for i in B)
        if A_measure and B_measure:
            boundary_measure = sum(v[i]*R[i,j] for i, j in product(A, B))
            A_connectivity = boundary_measure / A_measure
            B_connectivity = boundary_measure / B_measure
            connectivity = max(A_connectivity, B_connectivity)
            if I2 is None or connectivity < I2:
                I2 = connectivity
    return I2

def get_cheeger_bounds(R, v):
    """
    @param R: a reversible rate matrix
    @param v: stationary distribution
    @return: low, cheeger, high
    """
    cheeger = get_cheeger_constant(R, v)
    low = 0.5 * (cheeger**2) / -min(np.diag(R))
    high = 2 * cheeger
    return low, cheeger, high

def get_pi_mi_t2_approx(Q, v, t):
    """
    Second order taylor expansion for the contribution of each joint entry.
    @param Q: parent independent rate matrix
    @param v: stationary distribution
    @param t: amount of time
    """
    n = len(v)
    # get the randomization rate
    a = -np.trace(Q) / (n-1)
    return 0.5 * math.exp(-2 * a * t) * (n - 1)

def get_pi_mi_t2_diag_approx(Q, v, t):
    """
    Second order taylor expansion for only some contributions.
    Contributions of off-diagonal entries are computed exactly.
    @param Q: parent independent rate matrix
    @param v: stationary distribution
    @param t: amount of time
    """
    n = len(v)
    # get the randomization rate
    a = -np.trace(Q) / (n-1)
    x = math.exp(-a*t)
    h = ndot(v, 1-v)
    adjustment = (x*(1 - 0.5*x - math.log(1-x)) + math.log(1-x))*h
    return get_pi_mi_t2_approx(Q, v, t) + adjustment

def get_pi_mi(Q, v, t):
    """
    Try to get a closed form for the mutual information.
    This is for a parent independent process.
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

def get_pi_transition_matrix(Q, v, t):
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

def process(fs):
    nstates = fs.nstates
    np.set_printoptions(linewidth=200)
    t = fs.t
    ### sample a random time
    ##time_mu = 0.01
    ##t = random.expovariate(1 / time_mu)
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
    if fs.parent_indep:
        Q = np.outer(np.ones(nstates), v)
        Q -= np.diag(np.sum(Q, axis=1))
        pi_rescaling_factor = max(np.diag(R) / np.diag(Q))
        Q *= pi_rescaling_factor
        Z = msimpl.get_fast_meta_f81_autobarrier(Q)
    # Construct a child-independent process
    # with the same expected rate
    # as the sampled process
    if fs.child_indep:
        C = np.outer(1/v, np.ones(nstates))
        C -= np.diag(np.sum(C, axis=1))
        ci_rescaling_factor = np.max(R / C)
        #expected_rate = -ndot(np.diag(R), v)
        #ci_rescaling_factor = expected_rate / (nstates*(nstates-1))
        #ci_rescaling_factor = expected_rate / (nstates*nstates)
        C *= ci_rescaling_factor
        Q = C
    if fs.bipartitioned:
        Q = msimpl.get_fast_meta_f81_autobarrier(R)
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
    print >> out, 'eigenvalues of R:', scipy.linalg.eigvals(R)
    print >> out
    print >> out, 'relaxation rate of R:',
    print >> out, sorted(np.abs(scipy.linalg.eigvals(R)))[1]
    print >> out
    print >> out, 'expected rate of R:', mrate.Q_to_expected_rate(R)
    print >> out
    print >> out, 'cheeger bounds of R:', get_cheeger_bounds(R, v)
    print >> out
    print >> out, 'randomization rate of R:', get_randomization_rate(R, v)
    print >> out
    candidates = [get_randomization_candidate(R, v, i) for i in range(nstates)]
    if np.allclose(get_randomization_rate(R, v), candidates):
        print >> out, 'all candidates are equal to this rate'
    else:
        print >> out, 'not all candidates are equal to this rate'
    print >> out
    print >> out, 'simplified rate matrix Q:'
    print >> out, Q
    print >> out
    if fs.parent_indep:
        print >> out, 'parent independent rescaling factor:'
        print >> out, pi_rescaling_factor
        print >> out
    if fs.child_indep:
        print >> out, 'child independent rescaling factor:'
        print >> out, ci_rescaling_factor
        print >> out
    print >> out, 'eigenvalues of Q:', scipy.linalg.eigvals(Q)
    print >> out
    print >> out, 'relaxation rate of Q:',
    print >> out, sorted(np.abs(scipy.linalg.eigvals(Q)))[1]
    print >> out
    print >> out, 'expected rate of Q:', mrate.Q_to_expected_rate(Q)
    print >> out
    print >> out, 'cheeger bounds of Q:', get_cheeger_bounds(Q, v)
    print >> out
    print >> out, 'randomization rate of Q:', get_randomization_rate(Q, v)
    print >> out
    candidates = [get_randomization_candidate(Q, v, i) for i in range(nstates)]
    if np.allclose(get_randomization_rate(Q, v), candidates):
        print >> out, 'all candidates are equal to this rate'
    else:
        print >> out, 'warning: not all candidates are equal to this rate'
    print >> out
    print >> out, 'E(rate) of Q divided by logical entropy:',
    print >> out, mrate.Q_to_expected_rate(Q) / ndot(v, 1-v)
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
    print >> out, 'stationary distn logical entropy:', ndot(v, 1-v)
    print >> out
    # 
    P_by_hand = get_pi_transition_matrix(Q, v, t)
    print >> out, 'simplified-process transition matrix computed by hand:'
    print >> out, P_by_hand
    print >> out
    print >> out, 'simplified-process transition matrix computed by expm:'
    print >> out, scipy.linalg.expm(Q*t)
    print >> out
    #
    print >> out, 'simplified-process m.i. by hand:'
    print >> out, get_pi_mi(Q, v, t)
    print >> out
    print >> out, 'simplified-process m.i. by expm:'
    print >> out, ctmcmi.get_expected_ll_ratio(Q, t)
    print >> out
    #
    print >> out, 'original process m.i. by expm:'
    print >> out, ctmcmi.get_expected_ll_ratio(R, t)
    print >> out
    #
    print >> out, 'stationary distn Shannon entropy:'
    print >> out, -ndot(v, np.log(v))
    print >> out
    #
    if fs.parent_indep:
        print >> out, 'approximate simplified process m.i. 2nd order approx:'
        print >> out, get_pi_mi_t2_approx(Q, v, t)
        print >> out
        print >> out, 'approximate simplified process m.i. "better" approx:'
        print >> out, get_pi_mi_t2_diag_approx(Q, v, t)
        print >> out
        print >> out, '"f81-ization plus barrier" of pure f81-ization:'
        print >> out, Z
        print >> out
    #
    return out.getvalue().rstrip()

def get_response_content(fs):
    return process(fs) + '\n'

