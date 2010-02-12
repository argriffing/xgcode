"""
This module defines some distributions related to resequencing.

The defined distributions are for hidden Markov models.
For each distribution observations can be sampled,
and observations can be assigned a likelihood.
Observations are ordered sequences of four integers
corresponding to reads of A, C, G, and T.
"""

import unittest
import random
import math

import numpy as np
import scipy.stats
import scipy.maxentropy

import Util


def get_homozygous_distributions(randomization_rate):
    """
    Each distribution is over four states.
    @param randomization_rate: the probability that a read is randomized
    @return: a list of four distributions
    """
    distributions = []
    for favorite in range(4):
        d = [randomization_rate/4.0]*4
        d[favorite] = 1 - 3*randomization_rate/4.0
        distributions.append(d)
    return distributions

def get_heterozygous_distributions(randomization_rate):
    """
    Each distribution is over four states.
    @param randomization_rate: the probability that a read is randomized
    @return: a list of six distributions
    """
    distributions = []
    for first_index in range(4):
        for second_index in range(first_index):
            d = [randomization_rate/4.0]*4
            d[first_index] = .5 - randomization_rate/4.0
            d[second_index] = .5 - randomization_rate/4.0
            distributions.append(d)
    return distributions


class Mixture:
    """
    This allows sampling and likelihood calculations for a mixture model.
    This class can act as a HMM hidden state.
    """

    def __init__(self, states, distribution):
        """
        @param states: a sequence of hidden states
        @param distribution: the distribution of the hidden states
        """
        # do some validation
        if not len(states):
            raise ValueError('no states were specified')
        if len(states) != len(distribution):
            msg = 'the number of states should match the distribution length'
            raise ValueError(msg)
        if not np.allclose(sum(distribution), 1):
            raise ValueError('expected the distribution to sum to 1.0')
        if min(distribution) < 0:
            msg = 'expected the distribution to be a stochastic vector'
            raise ValueError(msg)
        # store the arguments, leaving out states with zero probability
        self.states = [state for state, d in zip(states, distribution) if d]
        self.distribution = [d for d in distribution if d]
        # precompute part of the likelihood
        self.log_distribution = [math.log(p) if p else float('-inf')
                for p in self.distribution]

    def get_posterior_distribution(self, observation):
        log_likelihoods = [state.get_log_likelihood(observation)
                for state in self.states]
        weighted_lls = [ll + log_p
                for ll, log_p in zip(log_likelihoods, self.log_distribution)]
        obs_ll = scipy.maxentropy.logsumexp(weighted_lls)
        return [math.exp(ll - obs_ll) for ll in weighted_lls]

    def sample_observation(self):
        """
        @return: an observation sampled from the distribution
        """
        state = self.states[Util.random_weighted_int(self.distribution)]
        return state.sample_observation()

    def get_likelihood(self, observation):
        return math.exp(self.get_log_likelihood(observation))

    def get_log_likelihood(self, observation):
        log_likelihoods = [state.get_log_likelihood(observation)
                for state in self.states]
        if all(ll==float('-inf') for ll in log_likelihoods):
            return float('-inf')
        weighted_log_likelihoods = [ll + log_p
                for ll, log_p in zip(log_likelihoods, self.log_distribution)]
        return scipy.maxentropy.logsumexp(weighted_log_likelihoods)


class UniformMixture:
    """
    This allows sampling and likelihood calculations for a mixture model.
    Each component of the mixture is equally likely.
    This class can act as a HMM hidden state.
    """

    def __init__(self, states):
        """
        @param states: a sequence of hidden states
        """
        self.states = states

    def sample_observation(self):
        return random.choice(self.states).sample_observation()

    def get_likelihood(self, observation):
        return math.exp(self.get_log_likelihood(observation))

    def get_log_likelihood(self, observation):
        log_likelihoods = [state.get_log_likelihood(observation) for state in self.states]
        if all(ll==float('-inf') for ll in log_likelihoods):
            return float('-inf')
        log_likelihood = scipy.maxentropy.logsumexp(log_likelihoods) - math.log(len(self.states))
        return log_likelihood


class SinglePatternState:
    """
    This is for when a single pattern is expected.
    Useful states are mixtures of these states.
    """

    def __init__(self, distribution, expected_coverage):
        """
        @param distribution: the expected nucleotide distribution
        @param expected_coverage: the read coverage at a position is poisson distributed with this expectation
        """
        self.distribution = distribution
        self.expected_coverage = expected_coverage

    def sample_observation(self):
        """
        @return: a sample of counts in each A, C, G, T state
        """
        n = np.random.poisson(self.expected_coverage)
        return np.random.multinomial(n, self.distribution)

    def get_log_likelihood(self, observation):
        if len(observation) != 4:
            raise ValueError('expected the observation to be a vector of four integers')
        n = sum(observation)
        accum = 0
        accum += Util.poisson_log_pmf(n, self.expected_coverage)
        accum += Util.multinomial_log_pmf(self.distribution, observation)
        return accum

    def get_likelihood(self, observation):
        return math.exp(self.get_log_likelihood(observation))


class FlatState:
    """
    This is supposed to be a somewhat flat distribution.
    Each of the four counts is sampled independently according to a geometric distribution.
    The distribution of the sum of counts is negative binomial.
    """

    def __init__(self, expected_coverage):
        """
        @param expected_coverage: the read coverage at a position is has this expectation
        """
        self.expected_coverage = expected_coverage

    def sample_observation(self):
        """
        @return: a sample of counts in each A, C, G, T state
        """
        mu = self.expected_coverage / 4.0
        pr = 1/(mu+1)
        return tuple(scipy.stats.geom.rvs(pr, loc=-1, size=4))

    def get_log_likelihood(self, observation):
        if len(observation) != 4:
            raise ValueError('expected the observation to be a vector of four integers')
        mu = self.expected_coverage / 4.0
        pr = 1/(mu+1)
        return sum(Util.geometric_log_pmf(obs, pr) for obs in observation)

    def get_likelihood(self, observation):
        return math.exp(self.get_log_likelihood(observation))


class Homozygous(UniformMixture):
    def __init__(self, randomization_rate, expected_coverage):
        distributions = get_homozygous_distributions(randomization_rate)
        states = [SinglePatternState(d, expected_coverage) for d in distributions]
        UniformMixture.__init__(self, states)

class Heterozygous(UniformMixture):
    def __init__(self, randomization_rate, expected_coverage):
        distributions = get_heterozygous_distributions(randomization_rate)
        states = [SinglePatternState(d, expected_coverage) for d in distributions]
        UniformMixture.__init__(self, states)

class Overcovered(UniformMixture):
    def __init__(self, randomization_rate, expected_coverage):
        distributions = get_homozygous_distributions(randomization_rate) + get_heterozygous_distributions(randomization_rate)
        states = [FlatState(expected_coverage)] + [SinglePatternState(d, expected_coverage) for d in distributions]
        UniformMixture.__init__(self, states)


def _gen_observations(n):
    """
    This is a helper function for testing.
    @param n: the total number of counts
    """
    for i in range(n+1):
        for j in range(n+1-i):
            for k in range(n+1-i-j):
                yield i, j, k, n-i-j-k


class TestReadCoverage(unittest.TestCase):

    def test_gen_observations(self):
        """
        Validate the sequences of possible observations.
        """
        # Assert that the right number of sequences is generated.
        # See OEIS A000292.
        for i in range(5):
            expected = ((i+1)*(i+2)*(i+3))/6
            observed = len(list(_gen_observations(i)))
            self.assertEqual(expected, observed)
        # Get a sequence.
        observations = list(_gen_observations(5))
        # Assert that each observation in the sequence is unique.
        self.assertEqual(len(observations), len(set(observations)))
        # Assert that each observation sums to the right number.
        for obs in observations:
            self.assertEqual(sum(obs), 5)

    def test_homozygous_sample(self):
        state = Homozygous(.2, 10)
        observation = state.sample_observation()
        self.assertEqual(len(observation), 4)

    def test_homozygous_likelihood(self):
        state = Homozygous(.2, 5)
        for n in range(10):
            p_observed = sum(state.get_likelihood(observation) for observation in _gen_observations(n))
            p_expected = scipy.stats.poisson.pmf(n, 5)
            self.assertTrue(np.allclose(p_observed, p_expected), (n, p_observed, p_expected))

    def test_heterozygous_sample(self):
        state = Heterozygous(.2, 10)
        observation = state.sample_observation()
        self.assertEqual(len(observation), 4)

    def test_heterozygous_likelihood(self):
        state = Heterozygous(.2, 5)
        for n in range(10):
            p_observed = sum(state.get_likelihood(observation) for observation in _gen_observations(n))
            p_expected = scipy.stats.poisson.pmf(n, 5)
            self.assertTrue(np.allclose(p_observed, p_expected), (n, p_observed, p_expected))

    def test_flat_sample(self):
        state = FlatState(10)
        observation = state.sample_observation()
        self.assertEqual(len(observation), 4)
        for x in observation:
            self.assertTrue(0 <= x)

    def test_flat_likelihood(self):
        expected_coverage = 5
        state = FlatState(expected_coverage)
        for n in range(10):
            p_observed = sum(state.get_likelihood(observation) for observation in _gen_observations(n))
            pr = 1.0 / (1.0 + expected_coverage / 4.0)
            p_expected = scipy.stats.nbinom.pmf(n, 4, pr)
            self.assertTrue(np.allclose(p_observed, p_expected), (n, p_observed, p_expected))

    def test_overcovered_sample(self):
        state = Overcovered(.2, 5)
        observation = state.sample_observation()
        self.assertEqual(len(observation), 4)
        for x in observation:
            self.assertTrue(0 <= x)

    def test_overcovered_likelihood(self):
        state = Overcovered(.2, 4)
        for n in range(8):
            p_observed = sum(state.get_likelihood(observation) for observation in _gen_observations(n))
            self.assertTrue(0 <= p_observed <= 1)

    def test_likelihood_ratios(self):
        """
        Assert that the distributions do what we want them to do.
        That is, assert that the objective classification of observations
        matches our subjective classification.
        """
        # define parameters
        good_coverage = 20
        bad_coverage = 60
        error_rate = .1
        # define hidden states
        homozygous = Homozygous(error_rate, good_coverage)
        heterozygous = Heterozygous(error_rate, good_coverage)
        overcovered = Overcovered(error_rate, bad_coverage)
        # define observations
        obs_a = [19, 1, 0, 0]
        obs_b = [10, 7, 0, 0]
        obs_c = [150, 12, 2, 1]
        # assert that the inference we want has the best log likelihood
        self.assertTrue(homozygous.get_log_likelihood(obs_a) > heterozygous.get_log_likelihood(obs_a))
        self.assertTrue(homozygous.get_log_likelihood(obs_a) > overcovered.get_log_likelihood(obs_a))
        self.assertTrue(heterozygous.get_log_likelihood(obs_b) > homozygous.get_log_likelihood(obs_b))
        self.assertTrue(heterozygous.get_log_likelihood(obs_b) > overcovered.get_log_likelihood(obs_b))
        self.assertTrue(overcovered.get_log_likelihood(obs_c) > homozygous.get_log_likelihood(obs_c))
        self.assertTrue(overcovered.get_log_likelihood(obs_c) > heterozygous.get_log_likelihood(obs_c))

    def test_large_observation_sizes(self):
        """
        Assert that the stats module can handle large input.
        """
        # define parameters
        good_coverage = 20
        bad_coverage = 60
        error_rate = .1
        # define hidden states
        homozygous = Homozygous(error_rate, good_coverage)
        heterozygous = Heterozygous(error_rate, good_coverage)
        overcovered = Overcovered(error_rate, bad_coverage)
        models = (homozygous, heterozygous, overcovered)
        # define the large observation
        obs = (0, 0, 605, 2)
        # attempt to get a log likelihood for each model
        for model in models:
            log_likelihood = model.get_log_likelihood(obs)
            self.failIf(math.isnan(log_likelihood))

    def test_myopic_model(self):
        """
        Test a mixture with a degenerately single-minded component.
        """
        # define a hidden state that is a mixture
        components = (FlatState(10), SinglePatternState((.9, .1, 0, 0), 10))
        model = UniformMixture(components)
        # define an observation that is completely incompatible with the second component
        obs = (0, 0, 10, 2)
        # assert that the log likelihood is reasonable
        log_likelihood = model.get_log_likelihood(obs)
        self.failIf(math.isnan(log_likelihood))

    def test_failed_mixture(self):
        """
        Test a mixture where no component can explain the data at all.
        """
        components = (
            SinglePatternState((0.9, 0.1, 0.0, 0.0), 10),
            SinglePatternState((0.8, 0.1, 0.1, 0.0), 10))
        model = UniformMixture(components)
        obs = (0, 0, 10, 2)
        log_likelihood = model.get_log_likelihood(obs)
        self.failIf(math.isnan(log_likelihood))

    def test_weighted_mixture_model_compatibility(self):
        """
        Create a weighted mixture that is the same as a uniform mixture.
        The likelihoods should be the same because it is the same model.
        """
        states = [FlatState(1), FlatState(10)]
        model_a = Mixture(states, [0.5, 0.5])
        model_b = UniformMixture(states)
        observation = (0,1,2,3)
        likelihood_a = model_a.get_likelihood(observation)
        likelihood_b = model_b.get_likelihood(observation)
        self.assertAlmostEqual(likelihood_a, likelihood_b)
        log_likelihood_a = model_a.get_log_likelihood(observation)
        log_likelihood_b = model_b.get_log_likelihood(observation)
        self.assertAlmostEqual(log_likelihood_a, log_likelihood_b)

    def test_weighted_mixture_model_inequality(self):
        """
        Create mixture models that differ only in their mixing proportions.
        The likelihoods should differ in predictable ways.
        """
        states = [FlatState(1), FlatState(10)]
        mixture_a = Mixture(states, [0.5, 0.5])
        mixture_b = Mixture(states, [0.4, 0.6])
        observation = (1,2,3,4)
        likelihood_a = mixture_a.get_likelihood(observation)
        likelihood_b = mixture_b.get_likelihood(observation)
        self.assertTrue(likelihood_a < likelihood_b)

if __name__ == '__main__':
    unittest.main()
