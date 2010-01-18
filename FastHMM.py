"""
Get the posterior distributions when some observations are missing, and use caching.

Observations should be hashable.
"""

import unittest
import math
import functools

import numpy as np

import TransitionMatrix
import DiscreteEndpoint
import MissingHMM
import HMM


class Model:

    def __init__(self, transition_object, hidden_state_objects, cache_size=0):
        """
        @param transition_object: returns a transition probability given two states and a distance
        @param hidden_state_objects: a conformant list of hidden state objects
        @param cache_size: the number of observations that are cached
        """
        nhidden = len(hidden_state_objects)
        self.transition_object = transition_object
        self.hidden_state_objects = hidden_state_objects
        self.initial_distribution = [transition_object.get_stationary_probability(i) for i in range(nhidden)]
        self.cache_size = cache_size
        self.cache = {}

    def get_likelihoods(self, obs):
        """
        This function uses memoization.
        @param obs: an emitted state
        @return: a tuple of likelihoods
        """
        likelihoods = self.cache.get(obs, None)
        if likelihoods:
            return likelihoods
        likelihoods = tuple(m.get_likelihood(obs) for m in self.hidden_state_objects)
        if len(self.cache) < self.cache_size:
            self.cache[obs] = likelihoods
        return likelihoods

    def scaled_forward_durbin(self, observations, distances):
        """
        For more information about the basic algorithm see HMM.scaled_forward_durbin.
        @param observations: the sequence of observations
        @param distances: positive integer distances between observations
        @return: the list of lists of scaled f variables, and the scaling variables
        """
        MissingHMM.validate_args(observations, distances)
        # precalculate
        nhidden = len(self.hidden_state_objects)
        nobs = len(observations)
        T = self.transition_object.get_transition_probability
        # initialize
        f = [[0]*nhidden for i in range(nobs)]
        s = [0]*nobs
        # define the initial unscaled f variable
        likelihoods = self.get_likelihoods(observations[0])
        for sink_index in range(nhidden):
            f[0][sink_index] = likelihoods[sink_index] * self.initial_distribution[sink_index]
        # define the positive scaling factor
        s[0] = sum(f[0])
        if not s[0]:
            raise ValueError('initial scaling factor is zero')
        # define the initial scaled f variable
        for sink_index in range(nhidden):
            f[0][sink_index] /= s[0]
        # define the subsequent f variables and scaling factors
        for i in range(1, nobs):
            obs = observations[i]
            likelihoods = self.get_likelihoods(obs)
            distance = distances[i-1]
            # define an unscaled f variable at this position
            for sink_index in range(nhidden):
                f[i][sink_index] = likelihoods[sink_index]
                p = 0
                for source_index in range(nhidden):
                    p += f[i-1][source_index] * T(source_index, sink_index, distance)
                f[i][sink_index] *= p
            # define the positive scaling factor at this position
            s[i] = sum(f[i])
            if not s[i]:
                raise ValueError('scaling factor is zero at position %d' % i)
            # define the scaled f variable at this position
            for sink_index in range(nhidden):
                f[i][sink_index] /= s[i]
        return f, s

    def scaled_backward_durbin(self, observations, distances, scaling_factors):
        """
        For more information about the basic algorithm see HMM.scaled_backward_durbin.
        The scaling factors must have been calculated using the scaled forward algorithm.
        @param observations: the sequence of observations
        @param distances: positive integer distances between observations
        @param scaling_factors: the scaling factor for each position
        @return: the list of lists of scaled b variables
        """
        MissingHMM.validate_args(observations, distances)
        # precalculate
        nhidden = len(self.hidden_state_objects)
        nobs = len(observations)
        T = self.transition_object.get_transition_probability
        # initialize
        b = [[0]*nhidden for i in observations]
        b[nobs-1] = [1/scaling_factors[nobs-1]]*nhidden
        for i in reversed(range(nobs-1)):
            obs = observations[i+1]
            likelihoods = self.get_likelihoods(obs)
            distance = distances[i]
            for source_index in range(nhidden):
                accum = 0
                for sink_index in range(nhidden):
                    p = 1.0
                    p *= T(source_index, sink_index, distance)
                    p *= likelihoods[sink_index]
                    p *= b[i+1][sink_index]
                    accum += p
                b[i][source_index] = accum / scaling_factors[i]
        return b

    def get_dp_info(self, observations, distances):
        """
        Do the dynamic programming and return the results.
        These results can be used for computing posterior expectations of
        emissions, of hidden states, and of transition counts.
        @param observations: the sequence of observations
        @param distances: positive integer distances between observations
        @return: observations, distances, f, s, b
        """
        MissingHMM.validate_args(observations, distances)
        f, s = self.scaled_forward_durbin(observations, distances)
        b = self.scaled_backward_durbin(observations, distances, s)
        return (observations, distances, f, s, b)

    def get_log_likelihood(self, dp_info):
        """
        @param dp_info: observations, distances, and forward, scaling, and backward variables
        @return: the log likelihood
        """
        # return sum of logs of scaling factors
        return sum(math.log(x) for x in dp_info[3])

    def scaled_posterior_durbin(self, dp_info):
        """
        @param dp_info: observations, distances, and forward, scaling, and backward variables
        @return: the list of position-specific posterior hidden state distributions
        """
        observations, distances, f, s, b = dp_info
        distributions = []
        for fs, bs, si in zip(f, b, s):
            distribution = [x*y*si for x, y in zip(fs, bs)]
            distributions.append(distribution)
        return distributions

    def scaled_ntransitions_expected(self, dp_info):
        """
        @param dp_info: observations, distances, and forward, scaling, and backward variables
        @return: the expected number of hidden state changes
        """
        observations, distances, f, s, b = dp_info
        nhidden = len(self.hidden_state_objects)
        nobs = len(observations)
        e_transitions = 0
        for i in range(1, nobs):
            obs = observations[i]
            distance = distances[i-1]
            likelihoods = self.get_likelihoods(obs)
            for sink_index in range(nhidden):
                for source_index in range(nhidden):
                    # compute the probability of a transition from the source to the sink
                    p = 1.0
                    p *= f[i-1][source_index]
                    p *= self.transition_object.get_transition_probability(source_index, sink_index, distance)
                    p *= likelihoods[sink_index]
                    p *= b[i][sink_index]
                    # add the number of expected transitions
                    e_transitions += p * self.transition_object.get_ntransitions_expected(source_index, sink_index, distance)
        return e_transitions


class DishonestCasino(Model):
    def __init__(self):
        fair_state = HMM.HiddenDieState(1/6.0)
        loaded_state = HMM.HiddenDieState(0.5)
        transition_matrix = np.array([[0.95, 0.05], [0.1, 0.9]])
        transition_object = TransitionMatrix.MatrixTransitionObject(transition_matrix)
        cache_size = 100
        Model.__init__(self, transition_object, [fair_state, loaded_state], cache_size)


class UniformState:
    def __init__(self, n):
        self.n = n
    def sample_observation(self):
        return random.randrange(self.n)
    def get_likelihood(self, observation):
        return 1.0 / self.n if observation in range(self.n) else 0.0
    def get_log_likelihood(self, observation):
        return -math.log(self.n) if observation in [0, 1] else float('-inf')


class FixedState:
    def __init__(self, k):
        self.k = k
    def sample_observation(self):
        return self.k
    def get_likelihood(self, observation):
        return 1.0 if observation == self.k else 0.0
    def get_log_likelihood(self, observation):
        return 0.0 if observation == self.k else float('-inf')


class TestFastHMM(unittest.TestCase):

    def test_scaled_posterior_durbin_compatibility(self):
        """
        Test the missing observation model when no observation is missing.
        """
        # define the models
        standard_hmm = HMM.DishonestCasino()
        missing_hmm = DishonestCasino()
        # define a sequence of observations
        observations = [1, 2, 6, 6, 1, 2, 3, 4, 5, 6]
        # define the (degenerate) distances between observations
        distances = [1]*(len(observations) - 1)
        # get posterior distributions with the standard algorithm
        standard_distributions = standard_hmm.scaled_posterior_durbin(observations)
        # get posterior distributions using the degenerate distances
        dp_info = missing_hmm.get_dp_info(observations, distances)
        missing_distributions = missing_hmm.scaled_posterior_durbin(dp_info)
        # assert that the distributions are the same
        self.assertTrue(np.allclose(standard_distributions, missing_distributions))

    def test_scaled_posterior_durbin_general(self):
        """
        This test is not strict but just checks some inequalities.
        """
        hmm = DishonestCasino()
        observations = [1, 2, 6, 6, 1, 2, 3, 4, 5, 6]
        distances_a = [1, 1, 1, 1, 1, 1, 1, 1, 1]
        distances_b = [1, 1, 1, 2, 1, 1, 1, 1, 1]
        # get posterior distributions
        distributions_a = hmm.scaled_posterior_durbin(hmm.get_dp_info(observations, distances_a))
        distributions_b = hmm.scaled_posterior_durbin(hmm.get_dp_info(observations, distances_b))
        # at offset 3 the first distribution should be more fair than the second distribution
        self.assertTrue(distributions_a[3][0] > distributions_b[3][0])
        # at offset 4 the first distribution should be less fair than the second distribution
        self.assertTrue(distributions_a[4][0] < distributions_b[4][0])

    def test_scaled_posterior_durbin_missing(self):
        """
        Test nontrivial distances against the uncached HMM implementation.
        """
        fast_hmm = DishonestCasino()
        slow_hmm = MissingHMM.DishonestCasino()
        observations = [1, 2, 6, 6, 1, 2, 3, 4, 5, 6]
        distances = [1, 1, 1, 2, 1, 1, 1, 1, 1]
        # get posterior distributions
        dp_info_a = fast_hmm.get_dp_info(observations, distances)
        distributions_a = fast_hmm.scaled_posterior_durbin(dp_info_a)
        distributions_b = slow_hmm.scaled_posterior_durbin(observations, distances)
        # assert that the posterior distributions are the same
        self.assertTrue(np.allclose(distributions_a, distributions_b))

    def test_scaled_ntransitions_expected_a(self):
        """
        Test the expected number of transitions given some partial observations.
        In this test there is no uncertainty about emitted states.
        """
        prandom = 0.1
        observations = [0, 0, 1]
        distances = [5, 7]
        states = [FixedState(i) for i in range(4)]
        nstates = len(states)
        # initialize the hidden markov model
        transition_object = TransitionMatrix.UniformTransitionObject(prandom, nstates)
        hmm = Model(transition_object, states)
        # compute the expected number of transitions using the hmm
        e_hmm = hmm.scaled_ntransitions_expected(hmm.get_dp_info(observations, distances))
        # compute the expected number of transitions directly
        e_direct = 0
        e_same, e_different = DiscreteEndpoint.get_expected_transitions_binomial(prandom, nstates, distances[0])
        e_direct += e_same
        e_same, e_different = DiscreteEndpoint.get_expected_transitions_binomial(prandom, nstates, distances[1])
        e_direct += e_different
        # assert that the results are the same
        self.assertAlmostEqual(e_hmm, e_direct)

    def test_scaled_ntransitions_expected_b(self):
        """
        Test the expected number of transitions given some partial observations.
        In this test there is complete randomization for hidden state transitions.
        """
        prandom = 1.0
        observations = [0, 0, 1]
        distances = [5, 7]
        states = [FixedState(i) for i in range(4)]
        nstates = len(states)
        # initialize the hidden markov model
        transition_object = TransitionMatrix.UniformTransitionObject(prandom, nstates)
        hmm = Model(transition_object, states)
        # compute the expected number of transitions using the hmm
        e_hmm = hmm.scaled_ntransitions_expected(hmm.get_dp_info(observations, distances))
        # compute the expected number of transitions directly
        e_direct = sum(distances) * (nstates - 1) / float(nstates)
        # assert that the results are the same
        self.assertAlmostEqual(e_hmm, e_direct)

    def test_scaled_ntransitions_expected_c(self):
        """
        Test the expected number of transitions given some partial observations.
        In this test there is complete uncertainty about emitted states.
        """
        prandom = 0.5
        observations = [0, 0, 1]
        distances = [5, 7]
        states = [UniformState(2) for i in range(4)]
        nstates = len(states)
        # initialize the hidden markov model
        transition_object = TransitionMatrix.UniformTransitionObject(prandom, nstates)
        hmm = Model(transition_object, states)
        # compute the expected number of transitions using the hmm
        e_hmm = hmm.scaled_ntransitions_expected(hmm.get_dp_info(observations, distances))
        # compute the expected number of transitions directly
        ptransition = prandom * (nstates - 1) / float(nstates)
        e_direct = sum(distances) * ptransition
        # assert that the results are the same
        self.assertAlmostEqual(e_hmm, e_direct)

    def test_scaled_ntransitions_expected_compatibility(self):
        fair_state = HMM.HiddenDieState(1/6.0)
        loaded_state = HMM.HiddenDieState(0.5)
        states = [fair_state, loaded_state]
        prandom = 0.1
        # define the old hmm
        transition_matrix = TransitionMatrix.get_uniform_transition_matrix(prandom, len(states))
        old_hmm = HMM.TrainedModel(transition_matrix, states)
        # define the new hmm
        cache_size = 100
        transition_object = TransitionMatrix.UniformTransitionObject(prandom, len(states))
        new_hmm = Model(transition_object, [fair_state, loaded_state], cache_size)
        # define a sequence of observations
        observations = [1, 2, 6, 6, 1, 2, 3, 4, 5, 6]
        # define the (degenerate) distances between observations
        distances = [1]*(len(observations) - 1)
        # use the old algorithm to get the expected number of transitions
        e_initial, A = old_hmm.scaled_transition_expectations_durbin(observations)
        ntransitions_expected_old = np.sum(A) - np.sum(np.diag(A))
        # use the new algorithm to get the expected number of transitions
        dp_info = new_hmm.get_dp_info(observations, distances)
        ntransitions_expected_new = new_hmm.scaled_ntransitions_expected(dp_info)
        # assert that the expected number of transitions are almost the same
        self.assertAlmostEqual(ntransitions_expected_old, ntransitions_expected_new)

    def test_scaled_ntransitions_expected_inequality(self):
        """
        In the dishonest casino a run of sixes means not much expected switching.
        Maybe this test should be moved to the HMM module.
        """
        fair_state = HMM.HiddenDieState(1/6.0)
        loaded_state = HMM.HiddenDieState(0.5)
        states = [fair_state, loaded_state]
        prandom = 0.1
        # define the new hmm
        cache_size = 100
        transition_object = TransitionMatrix.UniformTransitionObject(prandom, len(states))
        hmm = Model(transition_object, [fair_state, loaded_state], cache_size)
        # define sequences of observations
        observations_a = [6, 6, 6, 6, 1, 2, 2, 4, 5, 4]
        observations_b = [6, 6, 6, 6, 1, 6, 6, 6, 5, 6]
        # define the (degenerate) distances between observations
        distances = [1]*(len(observations_a) - 1)
        # use the algorithm to get the expected number of transitions for each observation sequence
        e_a = hmm.scaled_ntransitions_expected(hmm.get_dp_info(observations_a, distances))
        e_b = hmm.scaled_ntransitions_expected(hmm.get_dp_info(observations_b, distances))
        # assert that we see what we expect
        self.assertTrue(e_a > e_b)
        self.assertNotAlmostEqual(e_a, e_b)


if __name__ == '__main__':
    unittest.main()

