"""
This module is for reproducing tables from some papers.
"""

import unittest
from itertools import combinations
import math

import numpy as np
import scipy
from scipy import linalg
from scipy import integrate
from scipy import special

import StatsUtil
import MatrixUtil
import bernoulli
import wfengine
import wrightfisher
import Util

def params_to_mutation_fitness(N, params):
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
    fitness = 1.0 - np.array([g0, g1, g2, 0]) / float(N)
    return mutation, fitness

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

def get_transition_matrix_slow(N_diploid, k, mutation, fit):
    """
    Mutation probabilities are away from a fixed state.
    @param N_diploid: diploid population size
    @param k: number of alleles e.g. 4 for A,C,G,T
    @param mutation: k by k matrix of per-generation mutation probabilities
    @param fit: sequence of k fitness values
    @return: a transition matrix
    """
    N = N_diploid * 2
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
            #pi, pj = wrightfisher.genic_diallelic(fit[i], fit[j], h, N-h)
            #s = fit[i] - fit[j]
            s = 1 - fit[j] / fit[i]
            pi, pj = wrightfisher.genic_diallelic(1.0, 1.0 - s, h, N-h)
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

def get_transition_matrix(N_diploid, k, mutation, fit):
    """
    Mutation probabilities are away from a fixed state.
    @param N_diploid: diploid population size
    @param k: number of alleles e.g. 4 for A,C,G,T
    @param mutation: k by k matrix of per-generation mutation probabilities
    @param fit: sequence of k fitness values
    @return: a transition matrix
    """
    N = N_diploid * 2
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
        s = 1 - fit[j] / fit[i]
        pblock = np.exp(wfengine.create_genic_diallelic(N_diploid, s))
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

def get_scaled_fixation_probabilities(gammas):
    """
    The scaling factor is the same for each allele.
    It is something like 2N.
    @param gammas: scaled selections near zero, positive is better
    @return: a matrix of fixation probabilities
    """
    k = len(gammas)
    F = np.zeros((k, k))
    for i, gi in enumerate(gammas):
        for j, gj in enumerate(gammas):
            if i == j:
                continue
            F[i, j] = bernoulli.bgf(gj - gi)
    return F

def kimura_sojourn_helper(a, x):
    """
    Computes (exp(ax) - 1) / (exp(a) - 1) and accepts a=0.
    @param a: a scaled selection, can take any value
    @param x: a proportion between 0 and 1
    @return: a nonnegative value
    """
    if not a:
        return x
    else:
        return math.expm1(a*x) / math.expm1(a)

def _approx(p0, g, n0, n1, m0, m1):
    """
    This is a large population approximation.
    @param p0: proportion of allele 0 at the site
    @param g: difference of scaled selection
    @param n0: count of allele 0 in sample from child population
    @param n1: count of allele 1 in sample from child population
    @param m0: proportional to an expected number of substitutions
    @param m1: proportional to an expected number of substitutions
    @return: probability of the (n0, n1) selection given an n0+n1 size sample
    """
    p1 = 1 - p0
    # From Eq. (9) and (16) in McVean and Charlesworth 1999.
    # Note that a+b=1 ...
    a = kimura_sojourn_helper(g, p0)
    b = kimura_sojourn_helper(-g, p1)
    coeff = (m0*a + m1*b) / (p0 * p1)
    # get a binomial probability
    p = Util.choose(n0+n1, n0) * (p0 ** n0) * (p1 ** n1)
    # return the scaled probability
    return coeff * p

def diallelic_approximation(N_small, g, m0, m1):
    """
    This is a large population approximation.
    """
    hist = np.zeros(N_small+1)
    for n0 in range(1, N_small):
        n1 = N_small - n0
        hist[n0] = integrate.quad(
                _approx, 0, 1, args=(g, n0, n1, m0, m1))[0]
    return hist[1:-1] / np.sum(hist[1:-1])

def _approx_b(p0, g, n0, n1):
    """
    This is experimental.
    It should be related to the original similarly named function,
    but it is not as pedagogically useful.
    """
    p1 = 1 - p0
    coeff = kimura_sojourn_helper(g, p0) / (p0 * p1)
    p = Util.choose(n0+n1, n0) * (p0 ** n0) * (p1 ** n1)
    return coeff * p

def _approx_c(p0, g, n0m1, n1m1):
    """
    This is experimental.
    It should be related to the original similarly named function,
    but it is not as pedagogically useful.
    """
    p1 = 1 - p0
    coeff = kimura_sojourn_helper(g, p0)
    p = (p0 ** n0m1) * (p1 ** n1m1)
    return coeff * p

def diallelic_approximation_b(N_small, g, m0, m1):
    """
    This is experimental.
    It should be the same as the original similarly named function,
    but it is not as pedagogically useful.
    It uses the idea that if xxx is a mixture model of xxxa and xxxb
    then the compound distribution xxx-binomial is
    a mixture of xxxa-binomial and xxxb-binomial.
    """
    hist = np.zeros(N_small+1)
    for n0 in range(1, N_small):
        n1 = N_small - n0
        hist[n0] += m0 * integrate.quad(_approx_b, 0, 1, args=(g, n0, n1))[0]
        hist[n0] += m1 * integrate.quad(_approx_b, 0, 1, args=(-g, n1, n0))[0]
    return hist[1:-1] / np.sum(hist[1:-1])

def diallelic_approximation_b2(N_small, g, m0, m1):
    """
    This is experimental.
    It should be the same as the original similarly named function,
    but it is not as pedagogically useful.
    It uses the idea that if xxx is a mixture model of xxxa and xxxb
    then the compound distribution xxx-binomial is
    a mixture of xxxa-binomial and xxxb-binomial.
    """
    hist_a = np.zeros(N_small+1)
    hist_b = np.zeros(N_small+1)
    for n0 in range(1, N_small):
        n1 = N_small - n0
        hist_a[n0] += integrate.quad(_approx_b, 0, 1, args=(g, n0, n1))[0]
        hist_b[n0] += integrate.quad(_approx_b, 0, 1, args=(-g, n0, n1))[0]
    hist = m0 * hist_a + m1 * hist_b[::-1]
    return hist[1:-1] / np.sum(hist[1:-1])

def diallelic_approximation_c(N_small, g, m0, m1):
    """
    This is experimental.
    It should be the same as the original similarly named function,
    but it is not as pedagogically useful.
    It uses the idea that if xxx is a mixture model of xxxa and xxxb
    then the compound distribution xxx-binomial is
    a mixture of xxxa-binomial and xxxb-binomial.
    """
    hist = np.zeros(N_small - 1)
    for n0 in range(N_small - 1):
        n1 = N_small - 2 - n0
        c = Util.choose(n0+n1+2, n0+1)
        hist[n0] += m0 * c * integrate.quad(
                _approx_c, 0, 1, args=(g, n0, n1))[0]
        hist[n0] += m1 * c * integrate.quad(
                _approx_c, 0, 1, args=(-g, n1, n0))[0]
    return hist / np.sum(hist)

def diallelic_d_helper(n0, n1, g):
    if not g:
        return n0 / float(n0 + n1)
    else:
        return (special.hyp1f1(n0, n0 + n1, g) - 1) / math.expm1(g)

def diallelic_approximation_d(N_small, g, m0, m1):
    """
    This is experimental.
    The numerical integration should be replaced
    by a call to the confluent hypergeometric function hyp1f1.
    See also
    http://functions.wolfram.com/HypergeometricFunctions/
    Hypergeometric1F1/03/01/04/01/ .
    Also
    www.cs.unc.edu/Research/Image/MIDAG/p01/biostat/Digital_1.pdf .
    Also a gsl implementation gsl_sf_hyperg_1F1_int in hyperg_1F1.c
    specifically hyperg_1F1_ab_posint for positive integers a and b.
    Also
    http://mathworld.wolfram.com/
    ConfluentHypergeometricFunctionoftheFirstKind.html
    """
    hist = np.zeros(N_small + 1)
    for n0 in range(1, N_small):
        n1 = N_small - n0
        prefix = scipy.comb(n0+n1, n0) * special.beta(n0, n1)
        hist[n0] += m0 * prefix * diallelic_d_helper(n0, n1, g)
        hist[n0] += m1 * prefix * diallelic_d_helper(n1, n0, -g)
    return hist[1:-1] / np.sum(hist[1:-1])

def get_large_population_approximation(n, k, gammas, M):
    """
    @param n: sample size
    @param k: number of alleles e.g. 4 for A, C, G, T
    @param gammas: scaled selections near zero, positive is better
    @param M: mutation rate matrix
    @return: xxx
    """
    # Approximate the fixation probabilities.
    F = get_scaled_fixation_probabilities(gammas)
    # Compute the rate matrix as the hadamard product
    # of mutation rates and fixation probabilities,
    # and adjust the diagonal.
    Q = M*F
    Q -= np.diag(np.sum(Q, 1))
    # This is kind of a hack,
    # I should just get the stationary distribution directly from Q
    # without the expm.
    v = MatrixUtil.get_stationary_distribution(linalg.expm(Q))
    # Get sample allele frequencies associated with the
    # transitions between the fixed states of alleles 0 and 1.
    m0 = v[0] * M[0,1]
    m1 = v[1] * M[1,0]
    g = gammas[1] - gammas[0]
    return diallelic_approximation(n, g, m0, m1)

def get_test_mutation_fitness():
    mutation = np.array([
        [0.6, 0.2, 0.1, 0.1],
        [0.5, 0.1, 0.3, 0.1],
        [0.2, 0.1, 0.3, 0.4],
        [0.2, 0.3, 0.1, 0.4]])
    fitness = 1.0 - np.array([0.1, 0.2, 0.3, 0.4])
    return mutation, fitness

class TestTransitionMatrix(unittest.TestCase):
    def test_row_sums(self):
        N = 20
        k = 4
        mutation, fitness = get_test_mutation_fitness()
        P = get_transition_matrix(N, k, mutation, fitness)
        MatrixUtil.assert_transition_matrix(mutation)
        MatrixUtil.assert_transition_matrix(P)
    def test_stationary_distribution(self):
        N = 20
        k = 4
        mutation, fitness = get_test_mutation_fitness()
        P = get_transition_matrix(N, k, mutation, fitness)
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
        N_diploid = 10
        k = 4
        mutation, fitness = get_test_mutation_fitness()
        P_slow = get_transition_matrix_slow(N_diploid, k, mutation, fitness)
        P_fast = get_transition_matrix(N_diploid, k, mutation, fitness)
        self.assertTrue(np.allclose(P_slow, P_fast))

class TestLargePopulationApproximation(unittest.TestCase):
    def test_diallelic_approximation_equivalence_b(self):
        N_small = 10
        g = 1.5
        m0 = 2.0
        m1 = 3.0
        ha = diallelic_approximation(N_small, g, m0, m1)
        hb = diallelic_approximation_b(N_small, g, m0, m1)
        self.assertTrue(np.allclose(ha, hb))
    def test_diallelic_approximation_equivalence_b2(self):
        N_small = 10
        g = 1.5
        m0 = 2.0
        m1 = 3.0
        ha = diallelic_approximation(N_small, g, m0, m1)
        hb2 = diallelic_approximation_b2(N_small, g, m0, m1)
        self.assertTrue(np.allclose(ha, hb2))
    def test_diallelic_approximation_equivalence_c(self):
        N_small = 10
        g = 1.5
        m0 = 2.0
        m1 = 3.0
        ha = diallelic_approximation(N_small, g, m0, m1)
        hc = diallelic_approximation_c(N_small, g, m0, m1)
        self.assertTrue(np.allclose(ha, hc))
    def test_diallelic_approximation_equivalence_d(self):
        for g in (-1.5, 0, 1.5):
            N_small = 10
            m0 = 2.0
            m1 = 3.0
            ha = diallelic_approximation(N_small, g, m0, m1)
            hd = diallelic_approximation_d(N_small, g, m0, m1)
            self.assertTrue(np.allclose(ha, hd))

if __name__ == '__main__':
    unittest.main()

