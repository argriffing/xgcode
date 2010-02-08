"""
Condition a Markov chain on its endpoints.
The Markov chain is discrete and first order.

Find the expected number of transitions.
"""

import unittest
import itertools
import math

import numpy as np

import Util
import TransitionMatrix
import HMM
import iterutils


def get_expected_transitions_brute(prandom, nstates, nsteps):
    """
    This function is for transition matrices defined by their size and a single parameter.
    Use brute force to compute transition expectations.
    This function returns two values.
    The first value is the expected number of transitions
    when the endpoints are the same.
    The second value is the expected number of transitions
    when the endpoints are different.
    @param prandom: the probability of randomization at each step
    @param nstates: the number of states in the chain
    @param nsteps: one fewer than the length of the sequence
    @return: (expected_t_same, expected_t_different)
    """
    # handle corner cases
    if not nsteps:
        return 0.0, float('nan')
    if nsteps == 1:
        return 0.0, 1.0
    if not prandom:
        return 0.0, float('nan')
    # precalculate stuff
    p_notrans = prandom / nstates + (1 - prandom)
    p_particular_trans = prandom / nstates
    p_any_trans = p_particular_trans * (nstates - 1)
    # initialize probabilities
    total_p_different = 0
    total_p_same = 0
    # initialize expectations
    e_same = 0
    e_different = 0
    # define expectations
    for sequence in itertools.product(range(nstates), repeat=nsteps+1):
        # Calculate the probability of the sequence
        # and the number of transitions.
        ntransitions = 0
        p = 1.0 / nstates
        for a, b in iterutils.pairwise(sequence):
            if a == b:
                p *= p_notrans
            else:
                p *= p_particular_trans
                ntransitions += 1
        # add to the expectation
        if sequence[0] == sequence[-1]:
            total_p_same += p
            e_same += p * ntransitions
        else:
            total_p_different += p
            e_different += p * ntransitions
    e_same /= total_p_same
    e_different /= total_p_different
    return e_same, e_different

def get_expected_transitions_binomial(prandom, nstates, nsteps):
    """
    This function is for transition matrices defined by their size and a single parameter.
    Use binomial coefficients to compute transition expectations.
    @param prandom: the probability of randomization at each step
    @param nstates: the number of states in the chain
    @param nsteps: one fewer than the length of the sequence
    @return: (expected_t_same, expected_t_different)
    """
    # handle corner cases
    if not nsteps:
        return 0.0, float('nan')
    if nsteps == 1:
        return 0.0, 1.0
    if not prandom:
        return 0.0, float('nan')
    # precalculate stuff
    p_notrans = prandom / nstates + (1 - prandom)
    p_any_trans = 1.0 - p_notrans
    # precalculate expected probability of each endpoint pair state
    prandom_total = 1 - (1 - prandom)**nsteps
    p_notrans_total = prandom_total / nstates + (1 - prandom_total)
    # initialize expectations
    e_same = 0
    e_different = 0
    # define expectations
    for ntrans in range(nsteps+1):
        p_ntrans = math.exp(Util.binomial_log_pmf(ntrans, nsteps, p_any_trans))
        p_same = (1 - (1 - nstates)**(1 - ntrans))/nstates
        e_same += p_same * p_ntrans * ntrans
        e_different += (1 - p_same) * p_ntrans * ntrans
    e_same /= p_notrans_total
    e_different /= (1 - p_notrans_total)
    return e_same, e_different


class Chain:
    """
    This is an endpoint constrained Markov chain.
    The conditional transition expectation matrix can be computed
    using the forward and backward algorithms with scaling.
    """

    def __init__(self, transition_object):
        """
        @param transition_object: returns a transition probability given two states and a distance
        """
        self.nstates = transition_object.get_nstates()
        self.transition_object = transition_object
        self.initial_distribution = [transition_object.get_stationary_probability(i) for i in range(self.nstates)]

    def forward(self, initial_state, final_state, nsteps):
        """
        @param initial_state: the first state in the sequence
        @param final_state: the last state in the sequence
        @param nsteps: the number of transitions in the sequence
        @return: the list of lists of scaled f variables, and the scaling variables
        """
        T = self.transition_object.get_transition_probability
        # initialize
        f = [[0]*self.nstates for i in range(nsteps+1)]
        s = [0]*(nsteps+1)
        # define the initial f variable and scaling factor
        for state in range(self.nstates):
            f[0][state] = 1.0 if state == initial_state else 0.0
        s[0] = 1.0
        # define the subsequent f variables and scaling factors
        for i in range(1, nsteps+1):
            # define an unscaled f variable at this position
            for sink_index in range(self.nstates):
                if i < nsteps or sink_index == final_state:
                    f[i][sink_index] = 1.0
                else:
                    f[i][sink_index] = 0.0
                p = 0
                for source_index in range(self.nstates):
                    p += f[i-1][source_index] * T(source_index, sink_index)
                f[i][sink_index] *= p
            # define the positive scaling factor at this position
            s[i] = sum(f[i])
            if not s[i]:
                raise ValueError('scaling factor is zero at position %d' % i)
            # define the scaled f variable at this position
            for sink_index in range(self.nstates):
                f[i][sink_index] /= s[i]
        return f, s

    def backward(self, final_state, nsteps, scaling_factors):
        """
        The scaling factors must have been calculated using the scaled forward algorithm.
        @param final_state: the last state in the sequence
        @param nsteps: the number of transitions in the sequence
        @param scaling_factors: the scaling factor for each position
        @return: the list of lists of scaled b variables
        """
        T = self.transition_object.get_transition_probability
        b = [[0]*self.nstates for i in range(nsteps+1)]
        b[nsteps] = [1/scaling_factors[nsteps]]*self.nstates
        for i in reversed(range(nsteps)):
            for source_index in range(self.nstates):
                accum = 0
                for sink_index in range(self.nstates):
                    if i + 1 < nsteps or sink_index == final_state:
                        p = 1.0
                        p *= T(source_index, sink_index)
                        p *= b[i+1][sink_index]
                        accum += p
                b[i][source_index] = accum / scaling_factors[i]
        return b

    def get_dp_info(self, initial_state, final_state, nsteps):
        """
        Do the dynamic programming and return the results.
        These results can be used for computing posterior expectations of
        emissions, of hidden states, and of transition counts.
        @param initial_state: the first state in the sequence
        @param final_state: the last state in the sequence
        @param nsteps: the number of transitions in the sequence
        @return: initial_state, final_state, nsteps, f, s, b
        """
        f, s = self.forward(initial_state, final_state, nsteps)
        b = self.backward(final_state, nsteps, s)
        return (initial_state, final_state, nsteps, f, s, b)

    def get_transition_expectations(self, dp_info):
        """
        @param dp_info: this is from get_dp_info
        @return: a matrix of expected transition counts
        """
        initial_state, final_state, nsteps, f, s, b = dp_info
        T = self.transition_object.get_transition_probability
        # initialize the matrix of expected counts
        A = np.zeros((self.nstates, self.nstates))
        # get the expected counts for each transition
        for sink_index in range(self.nstates):
            for source_index in range(self.nstates):
                for i in range(1, nsteps+1):
                    if i < nsteps or sink_index == final_state:
                        p = 1.0
                        p *= f[i-1][source_index]
                        p *= T(source_index, sink_index)
                        p *= b[i][sink_index]
                        A[source_index, sink_index] += p
        return A

    def get_transition_expectations_brute(self, initial_state, final_state, nsteps):
        """
        @return: a matrix of expected transition counts
        """
        T = self.transition_object.get_transition_probability
        # initialize the matrix of expected counts
        A = np.zeros((self.nstates, self.nstates))
        # compute the probability of observing the final state conditional on the first state
        p_total = T(initial_state, final_state, nsteps)
        # iterate over all possible sequences of missing states
        for missing_sequence in itertools.product(range(self.nstates), repeat=nsteps-1):
            sequence = [initial_state] + list(missing_sequence) + [final_state]
            # get the probability of observing this continuation of the initial state
            p = 1.0
            for a, b in iterutils.pairwise(sequence):
                p *= T(a, b)
            # add the weighted transitions of each type
            for a, b in iterutils.pairwise(sequence):
                A[a, b] += p
        # divide by the total probability so that the conditioning is correct
        A /= p_total
        return A

    def get_expectations(self, dp_info):
        """
        @param dp_info: this is from get_dp_info
        @return: an expectation for each state
        """
        initial_state, final_state, nsteps, f, s, b = dp_info
        v = np.zeros(self.nstates)
        for fs, bs, si in zip(f, b, s)[1:-1]:
            for state, (x, y) in enumerate(zip(fs, bs)):
                v[state] += x*y*si
        return v

    def get_expectations_brute(self, initial_state, final_state, nsteps):
        """
        Get the number of times each state was expected to occur between the initial and final positions.
        @return: an expectation for each state
        """
        T = self.transition_object.get_transition_probability
        # initialize the vector of expected counts
        v = np.zeros(self.nstates)
        # compute the probability of observing the final state conditional on the first state
        p_total = T(initial_state, final_state, nsteps)
        # iterate over all possible sequences of missing states
        for missing_sequence in itertools.product(range(self.nstates), repeat=nsteps-1):
            sequence = [initial_state] + list(missing_sequence) + [final_state]
            # get the probability of observing this continuation of the initial state
            p = 1.0
            for a, b in iterutils.pairwise(sequence):
                p *= T(a, b)
            # add the weighted transitions of each type
            for state in missing_sequence:
                v[state] += p
        # divide by the total probability so that the conditioning is correct
        v /= p_total
        return v


class TestDiscreteEndpoint(unittest.TestCase):

    def test_brute(self):
        prandom = 0.001
        nstates = 4
        nsteps = 6
        e_same, e_different = get_expected_transitions_brute(prandom, nstates, nsteps)

    def test_binomial(self):
        for prandom in (0.01, 0.99, 1.0):
            for nstates in range(2, 6):
                for nsteps in range(1, 4):
                    expected = get_expected_transitions_brute(prandom, nstates, nsteps)
                    observed = get_expected_transitions_binomial(prandom, nstates, nsteps)
                    self.assertTrue(np.allclose(expected, observed))


    def test_chain_expected_transitions_brute_compatibility(self):
        """
        Test a reduced-parameter transition matrix for which a method already exists.
        """
        # the following parameters define the sequence distribution
        prandom = 0.1
        nstates = 3
        nsteps = 6
        # create the chain object
        transition_object = TransitionMatrix.UniformTransitionObject(prandom, nstates)
        chain = Chain(transition_object)
        # compare the expected number of changes
        A = chain.get_transition_expectations_brute(0, 0, nsteps)
        e_same_a = np.sum(A) - np.sum(np.diag(A))
        A = chain.get_transition_expectations_brute(0, 1, nsteps)
        e_different_a = np.sum(A) - np.sum(np.diag(A))
        e_same_b, e_different_b = get_expected_transitions_brute(prandom, nstates, nsteps)
        # the methods should give identical results
        self.assertAlmostEqual(e_same_a, e_same_b)
        self.assertAlmostEqual(e_different_a, e_different_b)

    def test_chain_expected_transitions(self):
        """
        Compare brute force results to dynamic programming results.
        """
        # define the sequence distribution
        T = np.array([
            [.1, .6, .3],
            [.1, .1, .8],
            [.8, .1, .1]])
        nstates = len(T)
        nsteps = 5
        # create the chain object
        transition_object = TransitionMatrix.MatrixTransitionObject(T)
        chain = Chain(transition_object)
        # compare the matrices defining the transition expectations
        for i_state, e_state in itertools.product(range(nstates), repeat=2):
            # compute the brute force result and the dynamic programming result
            A_brute = chain.get_transition_expectations_brute(i_state, e_state, nsteps)
            A_dynamic = chain.get_transition_expectations(chain.get_dp_info(i_state, e_state, nsteps))
            # assert that for each method the sum of the expectations is equal to the number of steps
            self.assertAlmostEqual(np.sum(A_brute), nsteps)
            self.assertAlmostEqual(np.sum(A_dynamic), nsteps)
            # assert that the methods give the same result
            self.assertTrue(np.allclose(A_brute, A_dynamic))

    def test_chain_expectations(self):
        """
        Compare brute force results to dynamic programming results.
        """
        # define the sequence distribution
        T = np.array([
            [.1, .6, .3],
            [.1, .1, .8],
            [.8, .1, .1]])
        nstates = len(T)
        nsteps = 5
        # create the chain object
        transition_object = TransitionMatrix.MatrixTransitionObject(T)
        chain = Chain(transition_object)
        # compare the matrices defining the transition expectations
        for i_state, e_state in itertools.product(range(nstates), repeat=2):
            # compute the brute force result and the dynamic programming result
            v_brute = chain.get_expectations_brute(i_state, e_state, nsteps)
            v_dynamic = chain.get_expectations(chain.get_dp_info(i_state, e_state, nsteps))
            # assert that for each method the sum of the expectations is equal to the number of steps minus one
            self.assertAlmostEqual(np.sum(v_brute), nsteps-1)
            self.assertAlmostEqual(np.sum(v_dynamic), nsteps-1)
            # assert that the methods give the same result
            self.assertTrue(np.allclose(v_brute, v_dynamic))


if __name__ == '__main__':
    unittest.main()

