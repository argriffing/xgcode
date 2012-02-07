"""
Compare algebraic connectivity of regular graphs vs. induced subgraphs.

The graphs are unweighted.
"""

from StringIO import StringIO
import math
import itertools

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
    form_objects = []
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_detailed_balance_error(Q):
    """
    @param Q: a rate matrix
    @return: a number that should be near zero if detailed balance is satisfied
    """
    p = mrate.R_to_distn(Q)
    errors = []
    nstates = len(Q)
    for i in range(nstates):
        for j in range(nstates):
            error = p[i] * Q[i, j] - p[j] * Q[j, i]
            errors.append(error)
    return max(abs(x) for x in errors)

def get_rate_matrix_summary(Q):
    out = StringIO()
    Q_v = mrate.R_to_distn(Q)
    Q_r = mrate.Q_to_expected_rate(Q)
    Q_t = mrate.R_to_relaxation_time(Q)
    print >> out, 'rate matrix:'
    print >> out, Q
    print >> out
    print >> out, 'this should be near zero for detailed balance:'
    print >> out, get_detailed_balance_error(Q)
    print >> out
    print >> out, 'computed stationary distribution:'
    print >> out, Q_v
    print >> out
    print >> out, 'expected rate:'
    print >> out, Q_r
    print >> out
    print >> out, 'relaxation time'
    print >> out, Q_t
    print >> out
    print >> out, '(expected rate) * (relaxation time):'
    print >> out, Q_r * Q_t
    print >> out
    print >> out
    return out.getvalue().rstrip()

def get_unweighted_cycle(nstates):
    """
    @return: a rate matrix
    """
    A = np.zeros((nstates, nstates))
    for i in range(nstates):
        a, b = i, (i+1) % nstates
        A[a,b] = 1
        A[b,a] = 1
    Q = A - np.diag(np.sum(A, axis=1))
    return Q

def get_unweighted_path(nstates):
    """
    @return: a rate matrix
    """
    A = np.zeros((nstates, nstates))
    for i in range(nstates-1):
        a, b = i, (i+1) % nstates
        A[a,b] = 1
        A[b,a] = 1
    Q = A - np.diag(np.sum(A, axis=1))
    return Q

def get_complete_graph(nstates):
    """
    @return: a rate matrix
    """
    A = np.ones((nstates, nstates), dtype=float)
    for i in range(nstates):
        A[i,i] = 0
    Q = A - np.diag(np.sum(A, axis=1))
    return Q

def get_hypercube(d):
    """
    The number of states is (2 ** d).
    @param d: the number of dimensions
    """
    return get_hamming_graph(2, d)

def get_hamming_graph(nresidues, nsites):
    """
    The number of states is (nresidues ** nsites).
    @return: a rate matrix
    """
    nstates = nresidues**nsites
    A = np.zeros((nstates, nstates))
    for alpha in itertools.product(range(nresidues), repeat=nsites):
        for beta in itertools.product(range(nresidues), repeat=nsites):
            alpha_index = sum(alpha[i]*(nresidues ** i) for i in range(nsites))
            beta_index = sum(beta[i]*(nresidues ** i) for i in range(nsites))
            hamming_dist = sum(1 for a, b in zip(alpha, beta) if a != b)
            if hamming_dist == 1:
                A[alpha_index, beta_index] = 1
    Q = A - np.diag(np.sum(A, axis=1))
    return Q


def get_response_content(fs):
    out = StringIO()
    np.set_printoptions(linewidth=200)
    # define the rate matrices
    Q_cycle = get_unweighted_cycle(10)
    Q_path = get_unweighted_path(3)
    # show the rate matrix summaries
    print >> out, 'unweighted 10-state cycle rate matrix'
    print >> out
    print >> out, get_rate_matrix_summary(Q_cycle)
    print >> out
    print >> out
    print >> out, 'unweighted 2-state path rate matrix'
    print >> out
    print >> out, get_rate_matrix_summary(Q_path)
    print >> out
    print >> out
    print >> out, 'unweighted 3-state complete graph rate matrix'
    print >> out
    print >> out, get_rate_matrix_summary(get_complete_graph(3))
    print >> out
    print >> out
    print >> out, 'unweighted 4-state complete graph rate matrix'
    print >> out
    print >> out, get_rate_matrix_summary(get_complete_graph(4))
    print >> out
    print >> out
    print >> out, 'unweighted 3-cube graph rate matrix'
    print >> out
    print >> out, get_rate_matrix_summary(get_hypercube(3))
    print >> out
    print >> out
    print >> out, 'unweighted 4-cube graph rate matrix'
    print >> out
    print >> out, get_rate_matrix_summary(get_hypercube(4))
    print >> out
    print >> out
    return out.getvalue()

