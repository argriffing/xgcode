"""
This module defines some more distributions related to resequencing.

The defined distributions are for hidden Markov models.
For each distribution observations can be sampled,
and observations can be assigned a likelihood.
Observations are ordered sequences of five integers
corresponding to reads of (A, C, G, T, gap).
"""

import unittest
import random
import math
import itertools

import numpy as np
import scipy.misc
import scipy.stats
from scipy.special import gammaln

import Util
import ReadCoverage


def get_homozygous_distributions(n, randomization_rate):
    """
    Each distribution is over n states.
    @param n: the number of states per distribution
    @param randomization_rate: the probability that a read is randomized
    @return: a list of n distributions
    """
    distributions = []
    for favorite in range(n):
        d = [randomization_rate/n]*n
        d[favorite] += 1 - randomization_rate
        distributions.append(d)
    return distributions

def get_heterozygous_distributions(n, randomization_rate):
    """
    Each distribution is over n states.
    @param n: the number of states per distribution
    @param randomization_rate: the probability that a read is randomized
    @return: a list of (n choose 2) distributions
    """
    distributions = []
    for first_index in range(n):
        for second_index in range(first_index):
            d = [randomization_rate/n]*n
            d[first_index] += (1 - randomization_rate)/2
            d[second_index] += (1 - randomization_rate)/2
            distributions.append(d)
    return distributions


class SinglePatternState:
    """
    This is for when a single pattern is expected.
    Useful states are mixtures of these states.
    """

    def __init__(self, distribution, expected_coverage):
        """
        The input distribution is like a discrete nucleotide distribution.
        @param distribution: the expected distribution
        @param expected_coverage: the read coverage at a position is poisson distributed with this expectation
        """
        self.distribution = distribution
        self.expected_coverage = expected_coverage
        # precalculate part of the log likelihood
        self.log_expected_coverage = math.log(self.expected_coverage)
        self.log_distribution = [math.log(p) for p in distribution]

    def sample_observation(self):
        """
        @return: a sample of counts in each state
        """
        n = np.random.poisson(self.expected_coverage)
        return np.random.multinomial(n, self.distribution)

    def get_log_likelihood(self, observation):
        """
        @param observation: like a tuple of nonnegative counts of nucleotides.
        """
        nstates = len(self.distribution)
        if len(observation) != nstates:
            raise ValueError('expected a vector of %d integers' % nstates)
        n = sum(observation)
        if not self.expected_coverage:
            if n:
                return float('-inf')
            else:
                return 0
        for d, c in zip(self.distribution, observation):
            if c and not d:
                return float('-inf')
        accum = 0
        accum += n * self.log_expected_coverage
        accum -= self.expected_coverage
        for log_p, obs in zip(self.log_distribution, observation):
            if obs > 1:
                accum -= gammaln(obs+1)
            if obs:
                accum += obs * log_p
        return accum

    def get_likelihood(self, observation):
        return math.exp(self.get_log_likelihood(observation))


class FlatState:
    """
    This is supposed to be a somewhat flat distribution.
    Each of the counts is sampled independently according to a geometric distribution.
    The distribution of the sum of counts is negative binomial.
    """

    def __init__(self, nstates, expected_coverage):
        """
        @param nstates: the number of different possible states each with a geometrically distributed count
        @param expected_coverage: the read coverage at a position has this expectation
        """
        # store the arguments
        self.nstates = nstates
        self.expected_coverage = expected_coverage
        # precalculate part of the log likelihood
        self.mu = self.expected_coverage / float(self.nstates)
        self.pr = 1/(self.mu+1)
        if self.pr != 0.0:
            self.log_pr = math.log(self.pr)
        if self.pr != 1.0:
            self.log_not_pr = math.log(1.0 - self.pr)

    def sample_observation(self):
        """
        @return: a sample of counts in each state
        """
        return tuple(scipy.stats.geom.rvs(self.pr, loc=-1, size=self.nstates))

    def get_log_likelihood(self, observation):
        if len(observation) != self.nstates:
            raise ValueError('expected a vector of %d integers' % self.nstates)
        if self.pr == 0.0:
            return float('-inf')
        if self.pr == 1.0:
            if any(observation):
                return float('-inf')
            else:
                return 0
        return sum(observation) * self.log_not_pr + self.nstates * self.log_pr

    def get_likelihood(self, observation):
        return math.exp(self.get_log_likelihood(observation))


class Homozygous:

    def __init__(self, nstates, randomization_rate, expected_coverage):
        """
        @param nstates: number of things you could be homozygous for
        @param randomization_rate: probability of a read giving a random result
        @param expected_coverage: expected sum of counts in an observation
        """
        # do validation
        if not (2 <= nstates):
            raise ValueError('the number of states should be at least two')
        if not (0.0 < randomization_rate <= 1.0):
            raise ValueError('the randomization rate should be a positive probability')
        if not (0 <= expected_coverage):
            raise ValueError('the expected coverage should be positive')
        # store arguments
        self.nstates = nstates
        self.randomization_rate = randomization_rate
        self.expected_coverage = expected_coverage
        # precalculate some stuff
        self.distributions = get_homozygous_distributions(nstates, randomization_rate)
        self.log_nstates = math.log(nstates)
        self.log_expected_coverage = math.log(expected_coverage)
        self.log_coeff = math.log(1 + self.nstates / self.randomization_rate - self.nstates)
        self.per_n = self.log_expected_coverage + math.log(randomization_rate) - self.log_nstates

    def sample_observation(self):
        """
        @return: a sample of counts in each state
        """
        distribution = random.choice(self.distributions)
        n = np.random.poisson(self.expected_coverage)
        return np.random.multinomial(n, distribution)

    def get_log_likelihood(self, observation):
        """
        @param observation: like a tuple of nonnegative counts of nucleotides.
        """
        if len(observation) != self.nstates:
            raise ValueError('expected a vector of %d integers' % self.nstates)
        n = sum(observation)
        sum_log_factorials = sum(gammaln(obs+1) if obs > 1 else 0 for obs in observation)
        accum = 0
        accum -= self.expected_coverage
        accum += n * self.per_n
        accum -= sum_log_factorials
        log_scaled_likelihoods = [obs*self.log_coeff for obs in observation]
        accum += scipy.misc.logsumexp(log_scaled_likelihoods)
        log_likelihood = accum - self.log_nstates
        if math.isnan(log_likelihood):
            raise ValueError('nan')
        return log_likelihood

    def get_likelihood(self, observation):
        return math.exp(self.get_log_likelihood(observation))


class Heterozygous:

    def __init__(self, nstates, randomization_rate, expected_coverage):
        """
        @param nstates: number of things you could be heterozygous for
        @param randomization_rate: probability of a read giving a random result
        @param expected_coverage: expected sum of counts in an observation
        """
        # do validation
        if not (2 <= nstates):
            raise ValueError('the number of states should be at least two')
        if not (0.0 < randomization_rate <= 1.0):
            raise ValueError('the randomization rate should be a positive probability')
        if not (0 <= expected_coverage):
            raise ValueError('the expected coverage should be positive')
        # store arguments
        self.nstates = nstates
        self.randomization_rate = randomization_rate
        self.expected_coverage = expected_coverage
        # precalculate some stuff
        self.distributions = get_heterozygous_distributions(nstates, randomization_rate)
        self.log_expected_coverage = math.log(expected_coverage)
        self.log_coeff = math.log(1 + self.nstates / (2.0*self.randomization_rate) - self.nstates/2.0)
        self.per_n = self.log_expected_coverage + math.log(randomization_rate) - math.log(nstates)
        self.log_nstates_choose_two = math.log((nstates*(nstates-1))/2)

    def sample_observation(self):
        """
        @return: a sample of counts in each state
        """
        distribution = random.choice(self.distributions)
        n = np.random.poisson(self.expected_coverage)
        return np.random.multinomial(n, distribution)

    def get_log_likelihood(self, observation):
        """
        @param observation: like a tuple of nonnegative counts of nucleotides.
        """
        if len(observation) != self.nstates:
            raise ValueError('expected a vector of %d integers' % self.nstates)
        n = sum(observation)
        sum_log_factorials = sum(gammaln(obs+1) if obs > 1 else 0 for obs in observation)
        accum = 0
        accum -= self.expected_coverage
        accum += n * self.per_n
        accum -= sum_log_factorials
        log_scaled_likelihoods = [(a+b)*self.log_coeff for a, b in itertools.combinations(observation, 2)]
        accum += scipy.misc.logsumexp(log_scaled_likelihoods)
        log_likelihood = accum - self.log_nstates_choose_two
        if math.isnan(log_likelihood):
            raise ValueError('nan')
        return log_likelihood

    def get_likelihood(self, observation):
        return math.exp(self.get_log_likelihood(observation))


class CachedSuperstate:
    """
    This is an abstract base for the 'good' and 'bad' superstates.
    """
    def __init__(self, cache_size=0):
        self.get_annotated_log_likelihood = Util.Cache(self._get_annotated_log_likelihood, cache_size)
    def _get_annotated_log_likelihood(self, observation):
        raise NotImplementedError()
    def sample_observation(self):
        raise NotImplementedError()
    def get_maximum_posterior(self, observation):
        log_likelihood, maximum_posterior = self.get_annotated_log_likelihood(observation)
        return maximum_posterior
    def get_log_likelihood(self, observation):
        log_likelihood, maximum_posterior = self.get_annotated_log_likelihood(observation)
        return log_likelihood
    def get_likelihood(self, observation):
        log_likelihood, maximum_posterior = self.get_annotated_log_likelihood(observation)
        return math.exp(log_likelihood)
    def get_annotated_likelihood(self, observation):
        log_likelihood, maximum_posterior = self.get_annotated_log_likelihood(observation)
        return math.exp(log_likelihood), maximum_posterior


class Good(CachedSuperstate):

    def __init__(self, randomization_rate, expected_coverage, cache_size=0):
        """
        The number of residues includes a gap residue.
        @param randomization_rate: the probability that a read is randomized
        @param expected_coverage: the read expected read coverage at a good position
        @param cache_size: the number of observations to remember
        """
        CachedSuperstate.__init__(self, cache_size)
        nresidues = 5
        mutation_rate = 1e-6
        state_hom = Homozygous(nresidues, randomization_rate, expected_coverage)
        state_het = Heterozygous(nresidues, randomization_rate, expected_coverage)
        self.states = [state_hom, state_het]
        self.distribution = [1 - mutation_rate, mutation_rate]
        self.state_names = ['hom', 'het']
        # precalculate stuff for efficiency
        self.get_posterior_distribution = Util.Cache(self._get_posterior_distribution, cache_size)
        self.log_distribution = [math.log(p) for p in self.distribution]

    def _get_annotated_log_likelihood(self, observation):
        """
        @param observation: an observation
        @return: the log likelihood and the MAP index
        """
        log_likelihoods = [m.get_log_likelihood(observation) for m in self.states]
        weighted_log_likelihoods = [ll + log_p for ll, log_p in zip(log_likelihoods, self.log_distribution)]
        ll, index = max((ll, i) for i, ll in enumerate(weighted_log_likelihoods))
        log_likelihood = scipy.misc.logsumexp(weighted_log_likelihoods)
        if math.isnan(log_likelihood):
            raise ValueError('nan')
        return log_likelihood, index

    def _get_posterior_distribution(self, observation):
        """
        @param observation: an observation
        @return: the posterior mixture distribution as a stochastic vector as a numpy array
        """
        log_likelihoods = [m.get_log_likelihood(observation) for m in self.states]
        weighted_lls = [ll + log_p for ll, log_p in zip(log_likelihoods, self.log_distribution)]
        weighted_lls = [x - max(weighted_lls) for x in weighted_lls]
        v = np.exp(weighted_lls)
        v /= np.sum(v)
        for x in v.tolist():
            if math.isnan(x):
                raise ValueError('nan: %s' % log_likelihoods)
        return v

    def sample_observation(self):
        """
        @return: an observation sampled from the distribution
        """
        return self.states[Util.random_weighted_int(self.distribution)].sample_observation()


class Bad(CachedSuperstate):

    def __init__(self, randomization_rate, expected_coverage, cache_size=0):
        """
        The number of residues includes a gap residue.
        @param randomization_rate: the probability that a read is randomized
        @param expected_coverage: the read expected read coverage at a good position
        @param cache_size: the number of observations to remember
        """
        CachedSuperstate.__init__(self, cache_size)
        nresidues = 5
        low_coverage = 2
        state_hom = Homozygous(nresidues, randomization_rate, expected_coverage)
        state_het = Heterozygous(nresidues, randomization_rate, expected_coverage)
        state_hom2x = Homozygous(nresidues, randomization_rate, expected_coverage*2)
        state_het2x = Heterozygous(nresidues, randomization_rate, expected_coverage*2)
        state_low = FlatState(nresidues, 2)
        state_weird = FlatState(nresidues, expected_coverage*2)
        state_high = FlatState(nresidues, expected_coverage*20)
        self.states = [state_hom, state_het, state_hom2x, state_het2x, state_low, state_weird, state_high]
        self.state_names = ['hom', 'het', 'hom2x', 'het2x', 'low', 'weird', 'high']
        # precalculate stuff for efficiency
        self.get_posterior_distribution = Util.Cache(self._get_posterior_distribution, cache_size)
        self.log_nstates = math.log(len(self.states))

    def _get_annotated_log_likelihood(self, observation):
        """
        @param observation: an observation
        @return: the log likelihood and the MAP index
        """
        log_likelihoods = [m.get_log_likelihood(observation) for m in self.states]
        ll, index = max((ll, i) for i, ll in enumerate(log_likelihoods))
        log_likelihood = scipy.misc.logsumexp(log_likelihoods) - self.log_nstates
        if math.isnan(log_likelihood):
            raise ValueError('nan')
        return log_likelihood, index

    def _get_posterior_distribution(self, observation):
        """
        @param observation: an observation
        @return: the posterior mixture distribution as a stochastic vector as a numpy array
        """
        log_likelihoods = [m.get_log_likelihood(observation) for m in self.states]
        lls = [x - max(log_likelihoods) for x in log_likelihoods]
        v = np.exp(lls)
        v /= np.sum(v)
        for x in v.tolist():
            if math.isnan(x):
                raise ValueError('nan: %s' % log_likelihoods)
        return v

    def sample_observation(self):
        return random.choice(self.states).sample_observation()


class TestReadCoverageGap(unittest.TestCase):

    def test_homozygous_distribution_compatibility(self):
        p = .1
        target_distributions = ReadCoverage.get_homozygous_distributions(p)
        query_distributions = get_homozygous_distributions(4, p)
        self.assertTrue(np.allclose(target_distributions, query_distributions))

    def test_heterozygous_distribution_compatibility(self):
        p = .1
        target_distributions = ReadCoverage.get_heterozygous_distributions(p)
        query_distributions = get_heterozygous_distributions(4, p)
        self.assertTrue(np.allclose(target_distributions, query_distributions))

    def test_flat_model_compatibility(self):
        target_state = ReadCoverage.FlatState(10)
        query_state = FlatState(4, 10)
        observation = (1, 2, 3, 4)
        target_likelihood = target_state.get_likelihood(observation)
        query_likelihood = query_state.get_likelihood(observation)
        self.assertAlmostEqual(query_likelihood, target_likelihood)

    def test_single_pattern_compatibility(self):
        d = (.1, .2, .3, .4)
        coverage = 10
        target_state = ReadCoverage.SinglePatternState(d, coverage)
        query_state = SinglePatternState(d, coverage)
        observation = (1, 2, 3, 4)
        target_likelihood = target_state.get_likelihood(observation)
        query_likelihood = query_state.get_likelihood(observation)
        self.assertAlmostEqual(query_likelihood, target_likelihood)

    def test_homozygous_compatibility(self):
        p = 0.1
        coverage = 10
        target_state = ReadCoverage.Homozygous(p, coverage)
        query_state = Homozygous(4, p, coverage)
        observation = (1, 2, 3, 4)
        target_likelihood = target_state.get_likelihood(observation)
        query_likelihood = query_state.get_likelihood(observation)
        self.assertAlmostEqual(query_likelihood, target_likelihood)

    def test_heterozygous_compatibility(self):
        p = 0.1
        coverage = 10
        target_state = ReadCoverage.Heterozygous(p, coverage)
        query_state = Heterozygous(4, p, coverage)
        observation = (1, 2, 3, 4)
        target_likelihood = target_state.get_likelihood(observation)
        query_likelihood = query_state.get_likelihood(observation)
        self.assertAlmostEqual(query_likelihood, target_likelihood)

    def test_superstate_likelhood(self):
        p = 0.1
        coverage = 10
        states = [Good(p, coverage), Bad(p, coverage)]
        for obs in ((0, 0, 0, 0, 0), (1, 0, 0, 0, 0), (1, 2, 3, 4, 5)):
            for state in states:
                p = state.get_likelihood(obs)

    def test_superstate_posterior(self):
        p = 0.1
        coverage = 10
        states = [Good(p, coverage), Bad(p, coverage)]
        for obs in ((0, 0, 0, 0, 0), (1, 0, 0, 0, 0), (1, 2, 3, 4, 5)):
            for state in states:
                p = state.get_posterior_distribution(obs)


if __name__ == '__main__':
    unittest.main()
