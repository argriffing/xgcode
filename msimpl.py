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

def get_fast_meta_f81(R, A):
    """
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
    Pa = np.sum(v[A])
    Pb = np.sum(v[B])
    flow = np.dot(v[A], np.sum(Rab, axis=1))
    Q = np.zeros((nstates, nstates))
    Q[np.ix_(A, A)] = get_fast_f81(Raa)
    Q[np.ix_(B, B)] = get_fast_f81(Rbb)
    Q[np.ix_(A, B)] = np.outer(np.ones(na), v[B]) * flow / (Pa * Pb)
    Q[np.ix_(B, A)] = np.outer(np.ones(nb), v[A]) * flow / (Pa * Pb)
    Q -= np.diag(np.sum(Q, axis=1))
    return Q

