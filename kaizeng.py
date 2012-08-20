"""
This module is for reproducing tables from some papers.
"""

import unittest
from itertools import combinations
import math

import numpy as np
from scipy import linalg

import StatsUtil
import MatrixUtil
import wfengine

def params_to_mutation_selection(N, params):
    """
    @param N: haploid population size
    @param params: parameters estimated by max likelihood in the 2011 paper
    """
    # define the hardcoded number of alleles
    k = 4
    # unpack the params
    theta, ka, kb, g0, g1, g2 = params
    # Expand the parameters into a higher dimensional
    # representation of mutation and selection.
    mutation = np.zeros((k, k))
    for i in range(k):
        for j in range(i+1,k):
            mutation[i,j] = theta / float(2*N)
    for i, j in ((0,1), (0,3), (1,3)):
        mutation[j,i] = ka * mutation[i,j]
    for i, j in ((0,2), (1,2), (2,3)):
        mutation[j,i] = kb * mutation[i,j]
    mutation += np.eye(k) - np.diag(np.sum(mutation, axis=1))
    selection = -np.array([g0, g1, g2, 0]) / float(N)
    return mutation, selection

def gen_states(N, k):
    """
    Yield states as count lists each of length k.
    The ith state corresponds to the fixation of the ith allele, when i<k.
    @param N: haploid population size
    @param k: number of alleles e.g. 4 for A,C,G,T
    """
    for i in range(k):
        state = [0]*k
        state[i] = N
        yield state
    for i, j in combinations(range(k), 2):
        for h in range(1, N):
            state = [0]*k
            state[i] = h
            state[j] = N-h
            yield state

def get_transition_matrix_slow(N, k, mutation, selection):
    """
    Mutation probabilities are away from a fixed state.
    @param N: haploid population size
    @param k: number of alleles e.g. 4 for A,C,G,T
    @param mutation: k by k matrix of per-generation mutation probabilities
    @param selection: sequence of k selection values
    @return: a transition matrix
    """
    fit = 1.0 + np.array(selection)
    states = [tuple(s) for s in gen_states(N,k)]
    nstates = len(states)
    s_to_i = dict((s, i) for i, s in enumerate(states))
    P = np.zeros((nstates, nstates))
    # Add rows corresponding to transitions from population states
    # for which an allele is currently fixed in the population.
    for i in range(k):
        P[i, i] = mutation[i, i]
        for j in range(k):
            if i == j:
                continue
            state = [0]*k
            state[i] = N-1
            state[j] = 1
            P[i, s_to_i[tuple(state)]] = mutation[i, j]
    # Add rows corresponding to transitions from polymorphic population states.
    for i, j in combinations(range(k), 2):
        for h in range(1, N):
            state = [0]*k
            state[i] = h
            state[j] = N-h
            index = s_to_i[tuple(state)]
            # Compute each child probability of having allele j.
            xi = h * fit[i]
            xj = (N-h) * fit[j]
            pi = xi / float(xi + xj)
            # Add entries corresponding to fixation of an allele.
            P[index, i] = math.exp(StatsUtil.binomial_log_pmf(N, N, pi))
            P[index, j] = math.exp(StatsUtil.binomial_log_pmf(0, N, pi))
            # Add entries corresponding to transitions to polymorphic states.
            for hsink in range(1, N):
                sink_state = [0]*k
                sink_state[i] = hsink
                sink_state[j] = N-hsink
                sink_index = s_to_i[tuple(sink_state)]
                logp = StatsUtil.binomial_log_pmf(hsink, N, pi)
                P[index, sink_index] = math.exp(logp)
    return P

def get_transition_matrix(N, k, mutation, selection):
    """
    Mutation probabilities are away from a fixed state.
    @param N: haploid population size
    @param k: number of alleles e.g. 4 for A,C,G,T
    @param mutation: k by k matrix of per-generation mutation probabilities
    @param selection: sequence of k selection values
    @return: a transition matrix
    """
    fit = 1.0 + np.array(selection)
    states = [tuple(s) for s in gen_states(N,k)]
    nstates = len(states)
    s_to_i = dict((s, i) for i, s in enumerate(states))
    P = np.zeros((nstates, nstates))
    # Add rows corresponding to transitions from population states
    # for which an allele is currently fixed in the population.
    for i in range(k):
        P[i, i] = mutation[i, i]
        for j in range(k):
            if i == j:
                continue
            state = [0]*k
            state[i] = N-1
            state[j] = 1
            P[i, s_to_i[tuple(state)]] = mutation[i, j]
    # Define transition matrices within a single diallelic subspace.
    for bi, (i, j) in enumerate(combinations(range(k), 2)):
        # Compute log probabilities.
        fitv = np.array([fit[i], fit[j]])
        log_distns = np.zeros((N+1, 2))
        for h in range(0, N+1):
            distn = np.array([h, (N-h)]) * fitv
            log_distns[h] = np.log(distn / np.sum(distn))
        # Compute the diallelic absorbing transition matrix.
        pblock = np.exp(wfengine.expand_multinomials(N, log_distns))
        ibegin = k + (N-1)*bi
        iend = ibegin + N - 1
        # The first index of pblock corresponds to fixation of j,
        # and the last index of pblock corresponds to fixation of i.
        # Incorporate fixation probabilities given various
        # nontrivial diallelic frequencies.
        P[ibegin:iend, i] = pblock[1:-1, -1]
        P[ibegin:iend, j] = pblock[1:-1, 0]
        # Incorporate transition probabilities among
        # nontrivial diallelic frequencies.
        # Note that the first and last row of pblock
        # are completely ignored.  This is intentional.
        P[ibegin:iend, ibegin:iend] = pblock[1:-1, 1:-1]
    return P

def get_test_mutation_selection():
    mutation = np.array([
        [0.6, 0.2, 0.1, 0.1],
        [0.5, 0.1, 0.3, 0.1],
        [0.2, 0.1, 0.3, 0.4],
        [0.2, 0.3, 0.1, 0.4]])
    selection = [0.1, 0.2, 0.3, 0.4]
    return mutation, selection

class TestTransitionMatrix(unittest.TestCase):
    def test_row_sums(self):
        N = 20
        k = 4
        mutation, selection = get_test_mutation_selection()
        P = get_transition_matrix(N, k, mutation, selection)
        MatrixUtil.assert_transition_matrix(mutation)
        MatrixUtil.assert_transition_matrix(P)
    def test_stationary_distribution(self):
        N = 20
        k = 4
        mutation, selection = get_test_mutation_selection()
        P = get_transition_matrix(N, k, mutation, selection)
        nstates = len(P)
        # one way to compute the stationary distribution
        w, vl = linalg.eig(P, left=True, right=False)
        v_eig = vl[:,0]
        v_eig = v_eig / np.sum(v_eig)
        # another way to compute the stationary distribution
        b = np.zeros(nstates)
        A = P.T - np.eye(nstates)
        A[0] = np.ones(nstates)
        b[0] = 1
        v_solve = linalg.solve(A, b)
        # Check that the two ways to get the stationary distribution
        # both give the same answer.
        self.assertTrue(np.allclose(v_eig, v_solve))
    def test_fast_slow_equivalence(self):
        N = 20
        k = 4
        mutation, selection = get_test_mutation_selection()
        P_slow = get_transition_matrix_slow(N, k, mutation, selection)
        P_fast = get_transition_matrix(N, k, mutation, selection)
        self.assertTrue(np.allclose(P_slow, P_fast))

if __name__ == '__main__':
    unittest.main()

