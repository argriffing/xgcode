"""
This module has functions that approximate rate matrices.

The idea is to approximate general reversible rate matrices
by simpler rate matrices whose properties are easier to analyze,
but which retain analytical connections to the original matrices.
The stuff about barriers is related to isoperimetry-type inequalities.
"""

from StringIO import StringIO
import random
import math
import itertools
from itertools import product

import numpy as np
import scipy
from scipy import linalg

import Form
import FormOut
import mrate
import ctmcmi
import iterutils
import combobreaker
import MatrixUtil
from MatrixUtil import ndot

#FIXME
def uniformize_split(R):
    """
    @param R: general reversible rate matrix
    @return: uniformization rate, transition matrix
    """
    n = len(R)
    v = mrate.R_to_distn(R)
    # For each row of R, determine the randomization rate
    # required to get the observed rate away.
    rates = [-R[i,i] / (1 - v[i]) for i in range(n)]
    # r is the uniformization rate
    r = max(rates)
    P = np.zeros((n,n))
    # Probability of no change is exp(- rate * t)
    for i in range(n):
        p_random = math.exp(-(r + R[i,i]))
        P[i] = R[i]
        P[i, i] += (1 - p_random)
    return r, P

#FIXME
def uniformize_join(r, P):
    """
    @param r: uniformization rate
    @param P: transition matrix
    @return: rate matrix
    """
    n = len(P)
    v = mrate.P_to_distn(P)
    Q = np.zeros((n,n))
    for i in range(n):
        p_random = None


def get_f81_ization_rate(Q, v):
    """
    @param Q: reversible rate matrix
    @param v: stationary distribution
    @return: randomization rate
    """
    nstates = len(v)
    return max(get_f81_ization_rate_candidate(Q, v, i) for i in range(nstates))

def get_f81_ization_rate_candidate(Q, v, i):
    """
    @param Q: reversible rate matrix
    @param v: stationary distribution
    @param i: a state
    @return: the randomization rate for the current state
    """
    return -Q[i,i] / (1 - v[i])

def get_barrier(R):
    """
    Return the subset of vertices on one side of a strong barrier.
    A strong barrier is one that minimizes the ratio of the
    flow across the barrier to the logical entropy of the partition.
    An alternative characterization is that this barrier
    minimizes the randomization rate of the corresponding
    Markov-ized 2-state process.
    @param R: general reversible rate matrix
    @return: a vertex subset
    """
    n = len(R)
    v = mrate.R_to_distn(R)
    best_subset = None
    best_ratio = None
    for A_tuple in iterutils.powerset(range(n)):
        # define the vertex set and its complement
        A = set(A_tuple)
        B = set(range(n)) - A
        Pa = sum(v[i] for i in A)
        Pb = sum(v[i] for i in B)
        if Pa and Pb:
            flow = sum(v[i]*R[i,j] for i, j in product(A, B))
            ratio = flow / (Pa * Pb)
            if (best_ratio is None) or (ratio < best_ratio):
                best_ratio = ratio
                best_subset = A
    print best_ratio
    return set(best_subset)

def get_fast_targeted_f81(R, v):
    """
    @param R: a general reversible rate matrix
    @param v: a target distribution that is not necessarily that of R
    @return: an F81 rate matrix with the target stationary distribution
    """
    nstates = len(R)
    if nstates == 1:
        return np.array([[0]])
    Q = np.outer(np.ones(nstates), v)
    Q -= np.diag(np.sum(Q, axis=1))
    Q *= max(np.diag(R) / np.diag(Q))
    return Q

def get_fast_f81(R):
    """
    A mutual information bound is conjectured.
    In particular the returned matrix is conjectured to have,
    at all positive times t, at most as much mutual information
    as the original matrix.
    @param R: a general reversible rate matrix
    @return: an F81 approximation with the same stationary distribution
    """
    nstates = len(R)
    if nstates == 1:
        return np.array([[0]])
    v = mrate.R_to_distn(R)
    Q = np.outer(np.ones(nstates), v)
    Q -= np.diag(np.sum(Q, axis=1))
    Q *= max(np.diag(R) / np.diag(Q))
    return Q

def get_fast_meta_f81_autobarrier(R):
    A = get_barrier(R)
    return get_fast_meta_f81(R, A)

def get_fast_two_state(R, A):
    n = len(R)
    v = mrate.R_to_distn(R)
    A = sorted(set(A))
    B = sorted(set(range(n)) - set(A))
    flow = sum(v[i] * R[i, j] for i, j in product(A, B))
    pa = sum(v[A])
    pb = sum(v[B])
    Q = np.array([
        [0, flow / pa],
        [flow / pb, 0]])
    Q -= np.diag(np.sum(Q, axis=1))
    return Q

def get_fast_two_state_autobarrier(R):
    A = get_barrier(R)
    return get_fast_two_state(R, A)

def get_fast_meta_f81_b(R, A):
    """
    The following conjecture does not work.
    A mutual information bound is conjectured.
    In particular the returned matrix is conjectured to have,
    at all positive times t, at most as much mutual information
    as the original matrix.
    @param R: a general reversible rate matrix
    @param A: a set of vertices on one side of a barrier
    @return: an approximation with two meta states
    """
    nstates = len(R)
    v = mrate.R_to_distn(R)
    A = sorted(set(A))
    B = sorted(set(range(nstates)) - set(A))
    na = len(A)
    nb = len(B)
    Raa = R[np.ix_(A, A)]
    Rbb = R[np.ix_(B, B)]
    Rab = R[np.ix_(A, B)]
    Rba = R[np.ix_(B, A)]
    Pa = np.sum(v[A])
    Pb = np.sum(v[B])
    flow_ab = np.dot(v[A], np.sum(Rab, axis=1))
    flow_ba = np.dot(v[B], np.sum(Rba, axis=1))
    if not np.allclose(flow_ab, flow_ba):
        raise ValueError('unequal R flow')
    # these lines are experimental to correct the stationary distribution
    #brunt_n_B = [sum(v[a]*R[a, b] for a in A) for b in B]
    #brunt_v_B = brunt_n_B / sum(brunt_n_B)
    #brunt_n_A = [sum(v[b]*R[b, a] for b in B) for a in A]
    #brunt_v_A = brunt_n_A / sum(brunt_n_A)
    # end experimental
    Q = np.zeros((nstates, nstates))
    Ra = Raa - np.sum(Raa, axis=1)
    Rb = Rbb - np.sum(Rbb, axis=1)
    Q[np.ix_(A, A)] = get_fast_targeted_f81(Ra, v[A] / sum(v[A]))
    Q[np.ix_(B, B)] = get_fast_targeted_f81(Rb, v[B] / sum(v[B]))
    Q[np.ix_(A, B)] = np.outer(np.ones(na), v[B]) * flow_ab / (Pa * Pb)
    Q[np.ix_(B, A)] = np.outer(np.ones(nb), v[A]) * flow_ab / (Pa * Pb)
    # begin experimental replacement
    #Q[np.ix_(A, B)] = np.outer(np.ones(na), brunt_v_B) * flow_ab / (Pa * sum(brunt_v_B))
    #Q[np.ix_(B, A)] = np.outer(np.ones(nb), brunt_v_A) * flow_ab / (sum(brunt_v_A) * Pb)
    Q -= np.diag(np.sum(Q, axis=1))
    # end experimental
    #
    # diagnostics
    qv = mrate.R_to_distn(Q)
    q_Pa = np.sum(qv[A])
    q_Pb = np.sum(qv[B])
    if not np.allclose(q_Pa, Pa):
        raise ValueError('error in meta stationary distribution')
    Qab = Q[np.ix_(A, B)]
    Qba = Q[np.ix_(B, A)]
    q_flow_ab = np.dot(qv[A], np.sum(Qab, axis=1))
    q_flow_ba = np.dot(qv[B], np.sum(Qba, axis=1))
    if not np.allclose(q_flow_ab, q_flow_ba):
        raise ValueError('unequal Q flow')
    return Q

def get_fast_meta_f81(R, A):
    """
    This also failed.
    This is a faster decaying function.
    Maybe the conjecture will work for it.
    A mutual information bound is conjectured.
    In particular the returned matrix is conjectured to have,
    at all positive times t, at most as much mutual information
    as the original matrix.
    @param R: a general reversible rate matrix
    @param A: a set of vertices on one side of a barrier
    @return: an approximation with two meta states
    """
    n = len(R)
    v = mrate.R_to_distn(R)
    f81_rate = get_f81_ization_rate(R, v)
    A = sorted(set(A))
    B = sorted(set(range(n)) - set(A))
    na = len(A)
    nb = len(B)
    Raa = R[np.ix_(A, A)]
    Rbb = R[np.ix_(B, B)]
    Rab = R[np.ix_(A, B)]
    Rba = R[np.ix_(B, A)]
    Pa = np.sum(v[A])
    Pb = np.sum(v[B])
    # Here the flow is unidirectional flow.
    # The unidirectional flow is the same quantity in both directions.
    flow_ab = np.dot(v[A], np.sum(Rab, axis=1))
    flow_ba = np.dot(v[B], np.sum(Rba, axis=1))
    flow = flow_ab
    if not np.allclose(flow_ab, flow_ba):
        raise ValueError('unequal R flow')
    # Get the conditional stationary distributions on each side
    # of the bipartition.
    va = v[A] / Pa
    vb = v[B] / Pb
    # Start making the rate matrix.
    Q = np.zeros((n, n))
    for a in A:
        # The sum of rates out of state a,
        # including a virtual flow that actually goes into a,
        # should be the f81_rate.
        # The rates into A should be proportional to the
        # stationary probabilities of states in A.
        # The rates into B should be proportional to the
        # stationary probabilities of states in B.
        # The ratio of rates into B to rates into A
        # should be the related to the flow
        # and to Pa and Pb.
        for tb in B:
            Q[a, tb] = v[tb] * flow / (Pa * Pb)
        for ta in A:
            Q[a, ta] = v[ta] * f81_rate / Pa
            #Q[a, ta] = (v[ta] / Pa) * (f81_rate - flow / (Pa * Pb))
    for b in B:
        for ta in A:
            Q[b, ta] = v[ta] * flow / (Pa * Pb)
        for tb in B:
            Q[b, tb] = v[tb] * f81_rate / Pb
            #Q[b, tb] = (v[tb] / Pb) * (f81_rate - flow / (Pa * Pb))
    Q -= np.diag(np.sum(Q, axis=1))
    # diagnostics to check the induced 2-state meta process
    qv = mrate.R_to_distn(Q)
    if not np.allclose(sum(qv[A]), sum(v[A])):
        raise ValueError('error in meta stationary distribution')
    Qab = Q[np.ix_(A, B)]
    Qba = Q[np.ix_(B, A)]
    q_flow_ab = np.dot(qv[A], np.sum(Qab, axis=1))
    q_flow_ba = np.dot(qv[B], np.sum(Qba, axis=1))
    if not np.allclose(q_flow_ab, q_flow_ba):
        raise ValueError('unequal Q flow')
    return Q

