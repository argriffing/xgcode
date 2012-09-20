"""
Check power to detect recessiveness. [UNFINISHED]

Assume diploids (required for recessiveness to be meaningful)
with genomic sites which are evolving
according to independent and identical processes.
Assume that the sites are diallelic.
The relevant evolutionary process parameters are
the fitness ratio between the two alleles,
the ratio of expected mutation rates from fixed population states
to polymorphic population states between the two alleles,
and a parameter that controls the amount of recessiveness.
Use simulation to look for the power to detect recessiveness
by computing the likelihood ratio statistic comparing a smaller null model
which assumes additive selection
to an alternative model that allows one more parameter
to define the recessiveness.
"""

from StringIO import StringIO
import time
import math

import numpy as np
from scipy import optimize
from scipy import linalg

import Form
import FormOut
import MatrixUtil
import StatsUtil
import kaizeng
import wfengine
import wrightfisher

def get_form():
    """
    @return: the body of a form
    """
    return [
            Form.Integer('N_diploid', 'diploid population size',
                100, low=1, high=500),
            Form.Integer('nsites', 'number of independent sites per genome',
                10000, low=1),
            Form.Integer('nalleles', 'number of sampled alleles per site',
                10, low=1),
            Form.Float('fit_A', 'fitness of allele A',
                1.01, low_exclusive=0),
            Form.Float('fit_a', 'fitness of allele a',
                1.00, low_exclusive=0),
            Form.Float('mut_A', 'expected mutation rate towards A from a',
                0.01, low_inclusive=0, high_inclusive=1),
            Form.Float('mut_a', 'expected mutation rate towards a from A',
                0.01, low_inclusive=0, high_inclusive=1),
            Form.Float('h', 'dominance of allele A (0.5 if additive)', 0.75),
            ]

def get_form_out():
    return FormOut.Report()

#XXX unused
def params_to_mutation_selection(N, params):
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

def params_to_distn(M_large, M_small, X):
    # unpack the parameters
    params = X.tolist()
    a, b, c = params
    # decode the parameters
    mut_ratio = math.exp(a)
    fit_ratio = math.exp(b)
    h = 0.5 * (math.tanh(c) + 1.0)
    # get the distribution implied by the parameters
    return get_sample_distn(M_large, M_small, mut_ratio, fit_ratio, h)

class G:
    def __init__(self, M_large, M_small, observed_counts):
        self.M_large = M_large
        self.M_small = M_small
        self.observed_counts = observed_counts
    def __call__(self, X):
        """
        @param X: three encoded parameters
        @return: negative log likelihood
        """
        # unpack the params into a finite distribution
        v = params_to_distn(self.M_large, self.M_small, X)
        # return the negative log likelihood
        return -StatsUtil.multinomial_log_pmf(v, self.observed_counts)

class G_additive:
    def __init__(self, M_large, M_small, observed_counts):
        self.M_large = M_large
        self.M_small = M_small
        self.observed_counts = observed_counts
    def __call__(self, X):
        """
        @return: negative log likelihood
        """
        # unpack the params into a finite distribution
        a, b = X.tolist()
        mut_ratio = math.exp(a)
        fit_ratio = math.exp(b)
        h = 0.5
        v = get_sample_distn(
                self.M_large, self.M_small,
                mut_ratio, fit_ratio, h)
        # return the negative log likelihood
        return -StatsUtil.multinomial_log_pmf(v, self.observed_counts)

class G_fit_only:
    def __init__(self, M_large, M_small, observed_counts):
        self.M_large = M_large
        self.M_small = M_small
        self.observed_counts = observed_counts
    def __call__(self, X):
        """
        @return: negative log likelihood
        """
        # unpack the params into a finite distribution
        a, = X.tolist()
        fit_ratio = math.exp(a)
        mut_ratio = 1.0
        h = 0.5
        v = get_sample_distn(
                self.M_large, self.M_small,
                mut_ratio, fit_ratio, h)
        # return the negative log likelihood
        return -StatsUtil.multinomial_log_pmf(v, self.observed_counts)

# XXX for testing
def get_two_allele_distribution(N_big, N_small, f0, f1, f_subsample):
    """
    Assumes small genic selection.
    Assumes small mutation.
    The mutational bias does not affect the distribution.
    @param N_big: total number of alleles in the population
    @param N_small: number of alleles sampled from the population
    @param f0: fitness of allele 0
    @param f1: fitness of allele 1
    @param f_subsample: subsampling function
    @return: distribution over all non-fixed population states
    """
    # construct a transition matrix
    nstates = N_big + 1
    P = np.zeros((nstates, nstates))
    for i in range(nstates):
        p0, p1 = wrightfisher.genic_diallelic(f0, f1, i, N_big - i)
        if i == 0:
            P[i, 1] = 1.0
        elif i == N_big:
            P[i, N_big - 1] = 1.0
        else:
            for j in range(nstates):
                logp = StatsUtil.binomial_log_pmf(j, N_big, p0)
                P[i, j] = math.exp(logp)
    # find the stationary distribution
    v = MatrixUtil.get_stationary_distribution(P)
    MatrixUtil.assert_distribution(v)
    if not np.allclose(v, np.dot(v, P)):
        raise ValueError('expected a left eigenvector with eigenvalue 1')
    # return the stationary distribution conditional on dimorphism
    print v
    distn = f_subsample(v, N_small)
    return distn[1:-1] / np.sum(distn[1:-1])

def get_polymorphic_diallelic_state_space(n):
    M = np.zeros((n-1, 2), dtype=int)
    for i in range(n-1):
        M[i, 0] = i+1
        M[i, 1] = n - M[i, 0]
    return M

def get_sample_distn(M_large, M_small, mut_ratio, fit_ratio, h):
    """
    @param M_large: an array of microstate counts
    @param M_small: an array of microstate counts
    @param mut_ratio: positive mutation rate ratio
    @param fit_ratio: positive fitness ratio
    @param h: recessiveness or dominance between 0 and 1
    @return: a distribution over polymorphic samples
    """
    # N is the population haploid count
    # n is the sample count
    N = sum(M_large[0])
    n = sum(M_small[0])
    nstates_large = len(M_large)
    nstates_small = len(M_small)
    if N-1 != nstates_large:
        raise Exception('internal error unexpected population state space')
    if N % 2:
        raise Exception('expected an even number of haploids in the population')
    if n-1 != nstates_small:
        raise Exception('internal error unexpected sample state space')
    # Get the value of s to construct the diallelic transition matrix.
    # It is positive when the lowercase allele is less fit.
    s = 1.0 - fit_ratio
    # Get the diallelic transition matrix.
    # The first state is allele A fixation
    # and the last state is allele a fixation.
    N_diploid = N // 2
    P = np.exp(wfengine.create_diallelic_recessive(N_diploid, s, h))
    MatrixUtil.assert_transition_matrix(P)
    # Get the expected number of visits
    # to polymorphic states, starting from (1, N-1) and from (N-1, 1),
    # assuming that fixed states are absorbing.
    I = np.eye(N - 1)
    Q = P[1:-1, 1:-1]
    B = np.zeros((N - 1, 2))
    B[0,0] = 1
    B[-1,-1] = 1
    X = linalg.solve((I - Q).T, B)
    # At this point X has two columns, each giving an array of expectations.
    # The first column gives expectations starting from a single A allele.
    # The second column gives expectations starting from a single a allele.
    # To get the expectations defined by the mixture,
    # take the weighted sum of the columns.
    e_large = X[:,0] * mut_ratio + X[:,1]
    # Compute a multivariate hypergeometric reduction
    # from the population state space expectations
    # to the sample state space expectations.
    e_small = wfengine.reduce_hypergeometric(M_large, M_small, e_large)
    # Compute the distribution.
    return e_small / e_small.sum()

def get_response_content(fs):
    np.set_printoptions(linewidth=200)
    out = StringIO()
    # extract user-supplied parameters
    N_diploid = fs.N_diploid
    nsites = fs.nsites
    nalleles = fs.nalleles
    fit_A = fs.fit_A
    fit_a = fs.fit_a
    mut_A = fs.mut_A
    mut_a = fs.mut_a
    h = fs.h
    #
    mut_ratio = mut_A / mut_a
    fit_ratio = fit_A / fit_a
    #
    N_hap = 2 * N_diploid
    M_large = get_polymorphic_diallelic_state_space(N_hap)
    M_small = get_polymorphic_diallelic_state_space(nalleles)
    #
    v_small = get_sample_distn(M_large, M_small, mut_ratio, fit_ratio, h)
    #
    print >> out, 'actual mut_ratio:', mut_ratio
    print >> out, 'actual fit_ratio:', fit_ratio
    print >> out, 'actual h:', h
    print >> out, 'implied finite polymorphic diallelic distribution:'
    print >> out, v_small
    print >> out
    #
    # sample from this distribution
    counts = np.random.multinomial(nsites, v_small)
    # try to estimate the parameters
    X0 = np.zeros(3)
    g = G(M_large, M_small, counts)
    Xopt = optimize.fmin(g, X0)
    #
    a, b, c = Xopt.tolist()
    mut_ratio_hat = math.exp(a)
    fit_ratio_hat = math.exp(b)
    h_hat = 0.5 * (math.tanh(c) + 1.0)
    v_small_hat = get_sample_distn(
            M_large, M_small,
            mut_ratio_hat, fit_ratio_hat, h_hat)
    #
    print >> out, 'estim. mut_ratio:', mut_ratio_hat
    print >> out, 'estim. fit_ratio:', fit_ratio_hat
    print >> out, 'estim. h:', h_hat
    print >> out, 'implied finite polymorphic diallelic distribution:'
    print >> out, v_small_hat
    print >> out
    #
    # constrain to additive selection
    X0 = np.zeros(2)
    g = G_additive(M_large, M_small, counts)
    Xopt = optimize.fmin(g, X0)
    a, b = Xopt.tolist()
    mut_ratio_hat = math.exp(a)
    fit_ratio_hat = math.exp(b)
    h_hat = 0.5
    v_small_hat = get_sample_distn(
            M_large, M_small,
            mut_ratio_hat, fit_ratio_hat, h_hat)
    print >> out, '-- inference assuming additive selection (h = 0.5) --'
    print >> out, 'estim. mut_ratio:', mut_ratio_hat
    print >> out, 'estim. fit_ratio:', fit_ratio_hat
    print >> out, 'estim. h:', h_hat
    print >> out, 'implied finite polymorphic diallelic distribution:'
    print >> out, v_small_hat
    print >> out
    #
    # constrain to additive selection and equal expected mutation
    X0 = np.zeros(1)
    g = G_fit_only(M_large, M_small, counts)
    Xopt = optimize.fmin(g, X0)
    a, = Xopt.tolist()
    mut_ratio_hat = 1.0
    fit_ratio_hat = math.exp(a)
    h_hat = 0.5
    v_small_hat = get_sample_distn(
            M_large, M_small,
            mut_ratio_hat, fit_ratio_hat, h_hat)
    print >> out, '-- inference assuming additive selection and equal mut --'
    print >> out, 'estim. mut_ratio:', mut_ratio_hat
    print >> out, 'estim. fit_ratio:', fit_ratio_hat
    print >> out, 'estim. h:', h_hat
    print >> out, 'implied finite polymorphic diallelic distribution:'
    print >> out, v_small_hat
    print >> out
    #
    return out.getvalue()


