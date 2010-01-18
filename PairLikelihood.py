"""Analyze or sample molecular sequence pairs according to a reversible continuous time Markov chain (CTMC).
"""

import unittest
import StringIO
import math
import random

import Util
import Newick
import Fasta
import RateMatrix
import MatrixUtil
import F84

class PairLikelihoodError(Exception):
    pass

def get_log_likelihood(t, sequence_pair, substitution_model):
    """
    @param t: the elapsed time, not the distance, between the sequence pairs
    @param sequence_pair: a pair of state sequences
    @param substitution_model: an object that represents the continuous time Markov chain
    @return: the log likelihood or None if there is no likelihood
    """
    # unpack the sequences
    sequence_a, sequence_b = sequence_pair
    # validate the input
    if len(sequence_pair) != 2:
        raise PairLikelihoodError('expected a pair of sequences')
    if len(sequence_a) != len(sequence_b):
        raise PairLikelihoodError('expected the two sequences to have the same length')
    # define the transition probabilities for the given model over the given time
    transition_dict = substitution_model.get_dictionary_transition_matrix(t)
    # define the list of ordered states
    ordered_states = substitution_model.states
    # define the stationary state distribution for the given model
    state_to_probability = dict(zip(ordered_states, substitution_model.get_stationary_distribution()))
    # compute the log likelihood
    log_likelihood = 0
    for state_a, state_b in zip(*sequence_pair):
        # calculate the probability factor due to this position in the aligned pair
        p = 1.0
        # consider the likelihood term due to the stationary probability of an observed state
        p *= state_to_probability[state_a]
        # consider the likelihood term due to the transition probability over the given amount of time
        p *= transition_dict[(state_a, state_b)]
        # if the probability is zero then do not take the logarithm of zero
        if not p:
            return None
        # add the contribution of this column to the log likelihood of the sequence pair
        log_likelihood += math.log(p)
    return log_likelihood

def simulate_sequence_pair(t, substitution_model, n_sites):
    """
    @param t: the elapsed time, not the distance, between the sequence pairs
    @param substitution_model: an object that represents the continuous time Markov chain
    @param n_sites: the length of each simulated sequence
    @return: a simulated pair of state sequences
    """
    # define the list of ordered states
    ordered_states = substitution_model.states
    # define the transition probabilities for the given model over the given time
    transition_dict = substitution_model.get_dictionary_transition_matrix(t)
    # define the stationary state distribution for the given model
    weight_state_pairs = zip(substitution_model.get_stationary_distribution(), ordered_states)
    # simulate the two sequences
    sequence_a = []
    sequence_b = []
    for i in range(n_sites):
        # sample the residue at this site on one sequence using the stationary distribution
        initial_state = Util.weighted_choice(weight_state_pairs)
        sequence_a.append(initial_state)
        # sample the residue at this site on the other sequence using the conditional distribution
        weight_target_pairs = [(transition_dict[(initial_state, target)], target) for target in ordered_states]
        target_state = Util.weighted_choice(weight_target_pairs)
        sequence_b.append(target_state)
    return (sequence_a, sequence_b)


class TestPairLikelihood(unittest.TestCase):

    def test_basic_log_likelihood(self):
        """Assert that the log likelihood function does not give an error for simple input."""
        kappa = 0.5
        nt_dist = (0.25, 0.25, 0.25, 0.25)
        model = F84.create_rate_matrix(kappa, nt_dist)
        sequence_pair = ('AAAA', 'ACGT')
        t = 2.0
        result = get_log_likelihood(t, sequence_pair, model)

    def test_unequal_log_likelihood(self):
        """Assert that the log likelihood function gives greater log likelihood to the more likely choice."""
        # define a sequence pair with many more transitions than transversions
        sequence_pair = ('AAAAAAAAAACCCCCCCCCC', 'ATGGGGGGGGCGTTTTTTTT')
        uniform_nt_dist = (0.25, 0.25, 0.25, 0.25)
        more_likely_model = F84.create_rate_matrix(1.0, uniform_nt_dist)
        less_likely_model = F84.create_rate_matrix(0.0, uniform_nt_dist)
        greater_log_likelihood = get_log_likelihood(1.0, sequence_pair, more_likely_model)
        smaller_log_likelihood = get_log_likelihood(1.0, sequence_pair, less_likely_model)
        self.failUnless(greater_log_likelihood > smaller_log_likelihood)

    def test_basic_simulation(self):
        """Assert that the sequence pair simulation does not give an error for simple input."""
        kappa = 0.5
        nt_dist = (0.25, 0.25, 0.25, 0.25)
        model = F84.create_rate_matrix(kappa, nt_dist)
        t = 2.0
        n_sites = 100
        result = simulate_sequence_pair(t, model, n_sites)
        self.assertEqual(len(result), 2)
        self.assertEqual(len(result[0]), 100)
        self.assertEqual(len(result[1]), 100)
        result_state_set = set(result[0]+result[1])
        self.failUnless(result_state_set <= set('ACGT'))



def main():
    pass

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('--test', action='store_true', dest='test', default=False)
    options, args = parser.parse_args()
    if options.test:
        suite = unittest.TestLoader().loadTestsFromTestCase(TestPairLikelihood)
        unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        main()

