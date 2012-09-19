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

#XXX unused
class G:
    def __init__(self, N, observed_counts):
        """
        @param N: haploid population size
        """
        self.N = N
        self.observed_counts = observed_counts
    def __call__(self, X):
        """
        @param X: six params defining mutation and selection
        @return: negative log likelihood
        """
        # define the hardcoded number of alleles
        k = 4
        # unpack the params
        params = X.tolist()
        theta, ka, kb, g0, g1, g2 = params
        if any(x < 0 for x in (theta, ka, kb)):
            return float('inf')
        mutation, fitnesses = kaizeng.params_to_mutation_fitness(
                self.N, params)
        # get the transition matrix
        P = kaizeng.get_transition_matrix(self.N, k, mutation, fitnesses)
        v = MatrixUtil.get_stationary_distribution(P)
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
    # Get the value of s to construct the diallelic transition matrix.
    # It is positive when the lowercase allele is less fit.
    s = 1.0 - fit_a / fit_A
    # Get the diallelic transition matrix.
    # The first state is allele A fixation
    # and the last state is allele a fixation.
    #P = np.exp(wfengine.create_diallelic_recessive(N_diploid, s, h))
    P = np.exp(wfengine.create_genic_diallelic(N_diploid, s))
    MatrixUtil.assert_transition_matrix(P)
    # Get the expected number of visits
    # to polymorphic states, starting from (1, N-1) and from (N-1, 1),
    # assuming that fixed states are absorbing.
    N_hap = 2 * N_diploid
    I = np.eye(N_hap - 1)
    Q = P[1:-1, 1:-1]
    B = np.zeros((N_hap - 1, 2))
    B[0,0] = 1
    B[-1,-1] = 1
    X = linalg.solve((I - Q).T, B)
    # At this point X has two columns, each giving an array of expectations.
    # The first column gives expectations starting from a single A allele.
    # The second column gives expectations starting from a single a allele.
    # To get the expectations defined by the mixture,
    # take the weighted sum of the columns.
    e_large = X[:,0]*mut_A + X[:,1]*mut_a
    #e_large = X[:,0] + X[:,1]
    # Define the population diallelic state space.
    M_large = np.zeros((N_hap-1, 2), dtype=int)
    for i in range(N_hap-1):
        M_large[i, 0] = i+1
        M_large[i, 1] = N_hap - M_large[i, 0]
    # Define the sample diallelic state space.
    M_small = np.zeros((nalleles-1, 2), dtype=int)
    for i in range(nalleles-1):
        M_small[i, 0] = i+1
        M_small[i, 1] = nalleles - M_small[i, 0]
    # Compute a multivariate hypergeometric reduction
    # from the population state space expectations
    # to the sample state space expectations.
    e_small = wfengine.reduce_hypergeometric(M_large, M_small, e_large)
    # Compute the distributions.
    v_large = e_large / e_large.sum()
    v_small = e_small / e_small.sum()
    # print stuff
    print >> out, v_large
    print >> out, v_small
    """
    print >> out, get_two_allele_distribution(
            N_hap, nalleles, 1.0, 1.0, 
            StatsUtil.subsample_pmf_without_replacement)
    """
    """
    for i in range(nsamples):
        counts = np.random.multinomial(nsites, v)
        X0 = np.array(params)
        g = G(N, counts)
        Xopt = optimize.fmin(g, X0)
        arr.append(Xopt)
    print >> out, np.array(arr)
    """
    return out.getvalue()

