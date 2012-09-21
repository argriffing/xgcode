"""
Check power to detect recessiveness in a two site diallelic compensatory system.

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
                100000, low=1),
            Form.Integer('nalleles', 'number of sampled alleles per site',
                10, low=1),
            Form.Float('mutation_rate', 'expected mutations per generation',
                0.01, low_exclusive=0, high_exclusive=1),
            Form.Float('Ns', 'selection coefficient Ns',
                1.00, low_exclusive=0),
            Form.Float('h', 'dominance h (0.5 if additive)', 0.5,
                low_exclusive=0, high_exclusive=1),
            ]

def get_form_out():
    return FormOut.Report()

# XXX old
def params_to_distn(M_large, M_small, X):
    # unpack the parameters
    params = X.tolist()
    a, b, c = params
    # decode the parameters
    mut_ratio = math.exp(a)
    fit_ratio = math.exp(b)
    h = StatsUtil.expit(c)
    # get the distribution implied by the parameters
    return get_sample_distn(M_large, M_small, mut_ratio, fit_ratio, h)

class G:
    def __init__(self, N, n, observed_counts):
        self.N = N
        self.n = n
        self.observed_counts = observed_counts
    def __call__(self, X):
        """
        @param X: three encoded parameters
        @return: negative log likelihood
        """
        # unpack the params into a finite distribution
        a, b, c = X.tolist()
        mutation_rate = StatsUtil.expit(a)
        fitness_ratio = math.exp(b)
        h = StatsUtil.expit(c)
        v = get_sample_distn(
                self.N, self.n,
                mutation_rate, fitness_ratio, h)
        # return the negative log likelihood
        return -StatsUtil.multinomial_log_pmf(v, self.observed_counts)

class G_additive:
    def __init__(self, N, n, observed_counts):
        self.N = N
        self.n = n
        self.observed_counts = observed_counts
    def __call__(self, X):
        """
        @return: negative log likelihood
        """
        # unpack the params into a finite distribution
        a, b = X.tolist()
        mutation_rate = StatsUtil.expit(a)
        fitness_ratio = math.exp(b)
        h = 0.5
        v = get_sample_distn(
                self.N, self.n,
                mutation_rate, fitness_ratio, h)
        # return the negative log likelihood
        return -StatsUtil.multinomial_log_pmf(v, self.observed_counts)

# XXX old
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

# XXX old
def get_polymorphic_diallelic_state_space(n):
    M = np.zeros((n-1, 2), dtype=int)
    for i in range(n-1):
        M[i, 0] = i+1
        M[i, 1] = n - M[i, 0]
    return M

def get_sample_distn(N, n, mutation_rate, fitness_ratio, h):
    """
    @param N: haploid pop size
    @param n: allele sample size
    @param mutation_rate: mutation rate
    @param fitness_ratio: fitness ratio
    @param h: dominance parameter 0.5 when additive
    @return: a distribution over n+1 mono-/di-morphic sample states
    """
    s = 1.0 - fitness_ratio
    P = np.exp(wfengine.create_diallelic_recessive(N // 2, s, h))
    MatrixUtil.assert_transition_matrix(P)
    # allow mutation out of the fixed states
    P[0, 0] = 1.0 - mutation_rate
    P[0, 1] = mutation_rate
    P[-1, -1] = 1.0 - mutation_rate
    P[-1, -2] = mutation_rate
    MatrixUtil.assert_transition_matrix(P)
    # get the population stationary distribution
    v_large = MatrixUtil.get_stationary_distribution(P)
    # get the allele distribution
    v_small = StatsUtil.subsample_pmf_without_replacement(v_large, n)
    return v_small

def get_response_content(fs):
    np.set_printoptions(linewidth=200)
    out = StringIO()
    # extract user-supplied parameters
    N_diploid = fs.N_diploid
    nsites = fs.nsites
    nalleles = fs.nalleles
    mutation_rate = fs.mutation_rate
    Ns = fs.Ns
    h = fs.h
    #
    N_hap = 2 * N_diploid
    N = N_hap
    n = nalleles
    s = fs.Ns / float(N)
    fitness_ratio = 1 - s
    v_small = get_sample_distn(N, n, mutation_rate, fitness_ratio, h)
    # sample from this distribution
    counts = np.random.multinomial(nsites, v_small)
    #
    negloglik = -StatsUtil.multinomial_log_pmf(
            v_small, counts)
    #
    print >> out, 'actual mutation rate:', mutation_rate
    print >> out, 'actual fitness ratio:', fitness_ratio
    print >> out, 'actual h:', h
    print >> out, 'implied finite polymorphic diallelic distribution:'
    print >> out, v_small
    print >> out, 'negative log likelihood:', negloglik
    print >> out
    #
    # try to estimate the parameters
    X0 = np.array([
        StatsUtil.logit(mutation_rate),
        math.log(fitness_ratio),
        StatsUtil.logit(h),
        ], dtype=float)
    g = G(N, n, counts)
    Xopt = optimize.fmin(g, X0)
    #
    a, b, c = Xopt.tolist()
    mutation_rate_hat = StatsUtil.expit(a)
    fitness_ratio_hat = math.exp(b)
    h_hat = StatsUtil.expit(c)
    v_small_hat = get_sample_distn(
            N, n,
            mutation_rate_hat, fitness_ratio_hat, h_hat)
    #
    negloglik_alt = g(Xopt)
    #
    print >> out, 'estim. mutation rate:', mutation_rate_hat
    print >> out, 'estim. fitness ratio:', fitness_ratio_hat
    print >> out, 'estim. h:', h_hat
    print >> out, 'implied finite polymorphic diallelic distribution:'
    print >> out, v_small_hat
    print >> out, 'negative log likelihood:', negloglik_alt
    print >> out
    #
    # constrain to additive selection
    X0 = np.array([
        StatsUtil.logit(mutation_rate),
        math.log(fitness_ratio),
        ], dtype=float)
    g = G_additive(N, n, counts)
    Xopt = optimize.fmin(g, X0)
    a, b = Xopt.tolist()
    mutation_rate_hat = StatsUtil.expit(a)
    fitness_ratio_hat = math.exp(b)
    h_hat = 0.5
    v_small_hat = get_sample_distn(
            N, n,
            mutation_rate_hat, fitness_ratio_hat, h_hat)
    #
    negloglik_null = g(Xopt)
    #
    print >> out, '-- inference assuming additive selection (h = 0.5) --'
    print >> out, 'estim. mutation rate:', mutation_rate_hat
    print >> out, 'estim. fitness ratio:', fitness_ratio_hat
    print >> out, 'estim. h:', h_hat
    print >> out, 'implied finite polymorphic diallelic distribution:'
    print >> out, v_small_hat
    print >> out, 'negative log likelihood:', negloglik_null
    print >> out
    #
    D = 2*(negloglik_null - negloglik_alt)
    print >> out, 'likelihood ratio test statistic:', D
    print >> out, 'chi squared 1-df 0.05 significance threshold:', 3.84
    print >> out
    #
    """
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
    """
    #
    return out.getvalue()


