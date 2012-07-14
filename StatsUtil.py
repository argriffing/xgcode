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
