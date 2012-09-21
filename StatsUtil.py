"""
Statistics utility functions.

Some attention is paid to edge cases and numerical stability.
"""

import unittest
from math import log
from math import exp

import numpy as np
import scipy.stats
from scipy.special import gammaln

import Util


def assert_probability(p):
    if p < 0:
        raise ValueError('%s (< 0) is not a probability' % p)
    elif p > 1:
        raise ValueError('%s (> 1) is not a probability' % p)

def binomial_log_pmf(observed_n, max_n, p_success):
    assert_probability(p_success)
    if p_success == 0.0:
        if observed_n:
            return float('-inf')
        else:
            return 0.0
    elif p_success == 1.0:
        if observed_n == max_n:
            return 0.0
        else:
            return float('-inf')
    accum = 0
    accum += gammaln(max_n + 1)
    accum -= gammaln(observed_n + 1)
    accum -= gammaln((max_n - observed_n) + 1)
    accum += observed_n * log(p_success)
    accum += (max_n - observed_n) * log(1.0 - p_success)
    return accum

def geometric_log_pmf(observed_n, pr):
    """
    @param observed_n: the number of completed events
    @param pr: the probability of quitting
    """
    assert_probability(pr)
    if pr == 0.0:
        return float('-inf')
    if pr == 1.0:
        if observed_n:
            return float('-inf')
        else:
            return log(pr)
    return observed_n * log(1.0 - pr) + log(pr)

def poisson_log_pmf(observed_n, expected_n):
    if not expected_n:
        if observed_n:
            return float('-inf')
        else:
            return 0.0
    accum = 0
    accum += observed_n * log(expected_n)
    accum -= expected_n
    accum -= gammaln(observed_n+1)
    return accum

def multinomial_log_pmf(distribution, counts):
    """
    This should be in scipy.stats but it isn't.
    @param distribution: the distribution over classes
    @param counts: the observed counts over classes
    """
    # check for a degeneracy
    for d, c in zip(distribution, counts):
        if c and not d:
            return float('-inf')
    n = sum(counts)
    # initialize the log probability mass
    accum = 0
    # add the contribution of n to the multinomial coefficient
    if n > 1:
        accum += gammaln(n+1)
    # add the contribution of the counts to the multinomial coefficient
    accum -= sum(gammaln(count+1) for count in counts if count > 1)
    # add the contribution of probabilities
    for p, count in zip(distribution, counts):
        if count:
            accum += count * log(p)
    return accum

def sample_distn_pmf_with_replacement(fdistn, n):
    """
    Pick a success probability at random then try n times.
    The returned distribution is over the number of successes.
    Each row of the finite distribution has two entries.
    The first entry is the probability of success and the
    second entry is the probability of that probability of success.
    @param fdistn: a finite distribution over success probabilities
    @param n: the total number of attempts
    @return: distribution over success counts
    """
    count_distn = np.zeros(n+1)
    for row in fdistn:
        psuccess = row[0]
        pfail = 1 - psuccess
        pstate = row[1]
        for nsuccess in range(n+1):
            nfail = n - nsuccess
            p = (psuccess ** nsuccess) * (pfail ** nfail)
            count_distn[i] += pstate * Util.choose(n, nsuccess) * p
    return count_distn

def subsample_pmf_with_replacement(distn, n):
    """
    Pick a number of successes in the population at random then pick a sample.
    The sample is picked with replacement.
    @param distn: distribution over the number of successes in the population
    @param n: the total number of attempts
    @return: distribution over success counts
    """
    N = len(distn) - 1
    count_distn = np.zeros(n+1)
    for Ns, pstate in enumerate(distn):
        ps = Ns / float(N)
        pf = 1 - ps
        for ns in range(n+1):
            nf = n - ns
            p = (ps ** ns) * (pf ** nf)
            count_distn[ns] += pstate * Util.choose(n, ns) * p
    return count_distn

def subsample_pmf_without_replacement(distn, n):
    """
    Pick a number of successes in the population at random then pick a sample.
    The sample is picked without replacement.
    @param distn: distribution over the number of successes in the population
    @param n: the total number of attempts
    @return: distribution over success counts
    """
    N = len(distn) - 1
    denominator = Util.choose(N, n)
    count_distn = np.zeros(n+1)
    for Ns, pstate in enumerate(distn):
        Nf = N - Ns
        for ns in range(n+1):
            nf = n - ns
            numerator = Util.choose(Ns, ns) * Util.choose(Nf, nf)
            count_distn[ns] += (pstate * numerator) / denominator
    return count_distn


def multinomial_log_pmf_vectorized(n, distn, counts):
    return gammaln(n+1) - np.sum(gammaln(counts+1)) + np.sum(counts*np.log(distn))


# TODO this is available as scipy.special.logit in scipy devel
def logit(p):
    """
    @param p: a probability
    @return: an unrestricted floating point number
    """
    return log(p) - log(1.0 - p)

# TODO this is available as scipy.special.expit in scipy devel
def expit(alpha):
    """
    @param alpha: an unrestricted floating point number
    @return: a probability
    """
    return exp(alpha) / (exp(alpha) + 1)


class TestStatsUtil(unittest.TestCase):

    def test_poisson_log_pmf(self):
        observed_n = 60
        expected_n = 20
        likelihood = scipy.stats.poisson.pmf(observed_n, expected_n)
        expected = log(likelihood)
        observed = poisson_log_pmf(observed_n, expected_n)
        self.assertAlmostEqual(expected, observed)

    def test_geometric_log_pmf(self):
        obs = 5
        pr = 0.1
        scipy_result = log(scipy.stats.geom.pmf(obs, pr, loc=-1))
        util_result = geometric_log_pmf(obs, pr)
        self.assertAlmostEqual(scipy_result, util_result)

    def test_binomial_log_pmf(self):
        """
        Test the binomial log pmf using an example from the internet.
        Suppose a die is tossed 5 times.
        What is the probability of getting exactly 2 fours?
        """
        observed_n = 2
        max_n = 5
        p_success = 1.0 / 6.0
        log_p = binomial_log_pmf(observed_n, max_n, p_success)
        p_computed = exp(log_p)
        p_expected = 0.160751028807
        self.assertTrue(np.allclose(p_expected, p_computed))


if __name__ == '__main__':
    unittest.main()
