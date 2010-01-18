"""
Get the position-specific posterior distribution of hidden states when some observations are missing.

The forward and backward algorithms are used according to
Biological Sequence Analysis by Durbin et al.
For numerical stability, scaling is used.
With these HMMs it is possible to calculate
the log likelihood of the observed sequence by efficiently summing over
all possible hidden states.
It is also possible to calculate the posterior distribution
of hidden states at each position for which an observation is available,
given the sequence of observed states.
Something that is not so straightforward with these HMMs is calculating
the expected transition counts,
so this component of the Baum-Welch algorithm is not implemented.
"""

import unittest

import numpy as np

import TransitionMatrix
import HMM


def validate_args(observations, distances):
    """
    @param observations: the sequence of observations
    @param distances: positive integer distances between observations
    """
    if len(observations) - len(distances) != 1:
        raise ValueError('the number of observations should be one more than the number of distances')
    for distance in distances:
        if int(distance) != distance:
            raise ValueError('each distance should be a positive integer')
        if distance < 1:
            raise ValueError('each distance should be a positive integer')


class MissingHMM:

    def __init__(self, transition_matrix, hidden_state_objects):
        """
        @param transition_matrix: a numpy array that is a transition matrix among hidden states
        @param hidden_state_objects: a conformant list of hidden state objects
        """
        self.transition_matrix = transition_matrix
        self.hidden_state_objects = hidden_state_objects
        self.stationary_distribution = TransitionMatrix.get_stationary_distribution(transition_matrix)
        self.initial_distribution = self.stationary_distribution

    def scaled_forward_durbin(self, observations, distances):
        """
        For more information about the basic algorithm see HMM.scaled_forward_durbin.
        @param observations: the sequence of observations
        @param distances: positive integer distances between observations
        @return: the list of lists of scaled f variables, and the scaling variables
        """
        validate_args(observations, distances)
        nhidden = len(self.hidden_state_objects)
        nobs = len(observations)
        f = [[0]*nhidden for i in range(nobs)]
        s = [0]*nobs
        # define the initial unscaled f variable
        for sink_index, sink_state in enumerate(self.hidden_state_objects):
            f[0][sink_index] = sink_state.get_likelihood(observations[0]) * self.initial_distribution[sink_index]
        # define the initial scaling factor
        s[0] = sum(f[0])
        # define the initial scaled f variable
        for sink_index in range(nhidden):
            f[0][sink_index] /= s[0]
        # define the subsequent f variables and scaling factors
        for i in range(1, nobs):
            obs = observations[i]
            # define the position specific transition matrix
            T = np.linalg.matrix_power(self.transition_matrix, distances[i-1])
            # define an unscaled f variable at this position
            for sink_index, sink_state in enumerate(self.hidden_state_objects):
                f[i][sink_index] = sink_state.get_likelihood(obs)
                p = 0
                for source_index, source_state in enumerate(self.hidden_state_objects):
                    p += f[i-1][source_index] * T[source_index, sink_index]
                f[i][sink_index] *= p
            # define the scaling factor at this position
            s[i] = sum(f[i])
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
        validate_args(observations, distances)
        nhidden = len(self.hidden_state_objects)
        nobs = len(observations)
        b = [[0]*nhidden for i in observations]
        b[nobs-1] = [1/scaling_factors[nobs-1]]*nhidden
        for i in reversed(range(nobs-1)):
            T = np.linalg.matrix_power(self.transition_matrix, distances[i])
            for source_index, source_state in enumerate(self.hidden_state_objects):
                accum = 0
                for sink_index, sink_state in enumerate(self.hidden_state_objects):
                    p = 1.0
                    p *= T[source_index, sink_index]
                    p *= sink_state.get_likelihood(observations[i+1])
                    p *= b[i+1][sink_index]
                    accum += p
                b[i][source_index] = accum / scaling_factors[i]
        return b

    def scaled_posterior_durbin(self, observations, distances):
        """
        For more information about the basic algorithm see HMM.scaled_backward_durbin.
        @param observations: the sequence of observations
        @param distances: positive integer distances between observations
        @return: the list of position-specific posterior hidden state distributions
        """
        validate_args(observations, distances)
        f, s = self.scaled_forward_durbin(observations, distances)
        b = self.scaled_backward_durbin(observations, distances, s)
        distributions = []
        for fs, bs, si in zip(f, b, s):
            distribution = [x*y*si for x, y in zip(fs, bs)]
            distributions.append(distribution)
        return distributions


class DishonestCasino(MissingHMM):
    def __init__(self):
        fair_state = HMM.HiddenDieState(1/6.0)
        loaded_state = HMM.HiddenDieState(0.5)
        transition_matrix = np.array([[0.95, 0.05], [0.1, 0.9]])
        MissingHMM.__init__(self, transition_matrix, [fair_state, loaded_state])


class TestMissingHMM(unittest.TestCase):

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
        missing_distributions = missing_hmm.scaled_posterior_durbin(observations, distances)
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
        distributions_a = hmm.scaled_posterior_durbin(observations, distances_a)
        distributions_b = hmm.scaled_posterior_durbin(observations, distances_b)
        # at offset 3 the first distribution should be more fair than the second distribution
        self.assertTrue(distributions_a[3][0] > distributions_b[3][0])
        # at offset 4 the first distribution should be less fair than the second distribution
        self.assertTrue(distributions_a[4][0] < distributions_b[4][0])


if __name__ == '__main__':
    unittest.main()

