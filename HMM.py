"""
Do some Markov model stuff.
For now, use hidden states from a finite alphabet,
but allow emissions to be more complicated objects.
For example, allow the observed sequence
to be a sequence of fixed length vectors of unbounded integers.
"""

import random
import unittest
import math
import itertools

import numpy as np
import scipy.maxentropy

import Util
import TransitionMatrix


def get_example_rolls():
    """
    This is a helper function for testing.
    See figure 3.5 in the eighth printing of
    Biological Sequence Analysis.
    Lines alternate between rolls and estimates.
    Each roll is a die roll in (1, 2, 3, 4, 5, 6).
    Each estimate is a 'Fair' vs. 'Loaded' estimate.
    @return: (300 observations, 300 viterbi estimates)
    """
    lines = [
            '315116246446644245311321631164',
            'FFFFFFFFFFFFFFFFFFFFFFFFFFFFFF',
            '152133625144543631656626566666',
            'FFFFFFFFFFFFFFFFFFLLLLLLLLLLLL',
            '651166453132651245636664631636',
            'LLLLLLFFFFFFFFFFFFLLLLLLLLLLLL',
            '663162326455236266666625151631',
            'LLLLLLLLLLLLLLLLLLLLLLFFFFFFFF',
            '222555441666566563564324364131',
            'FFFFFFFFFFFFFFFFFFFFFFFFFFFFFF',
            '513465146353411126414626253356',
            'FFFFFFFFFFFFFFFFFFFFFFFFFFFFFL',
            '366163666466232534413661661163',
            'LLLLLLLLLLLLFFFFFFFFFFFFFFFFFF',
            '252562462255265252266435353336',
            'FFFFFFFFFFFFFFFFFFFFFFFFFFFFFF',
            '233121625364414432335163243633',
            'FFFFFFFFFFFFFFFFFFFFFFFFFFFFFF',
            '665562466662632666612355245242',
            'LLLLLLLLLLLLLLLLLLLFFFFFFFFFFF']
    observation_lines = [line for i, line in enumerate(lines) if i%2 == 0]
    estimate_lines = [line for i, line in enumerate(lines) if i%2 == 1]
    observations = [int(x) for x in ''.join(observation_lines)]
    estimates = [x for x in ''.join(estimate_lines)]
    return observations, estimates

def get_distribution(models, obs):
    """
    This assumes that the models are equally believed before any observation.
    This function does not use any Markov stuff.
    @param models: a list of sources of log likelihoods
    @param obs: an observation
    @return: the stochastic vector of relative posterior beliefs of the models as a numpy array
    """
    v = np.array([m.get_log_likelihood(obs) for m in models])
    v = np.exp(v - np.max(v))
    v /= np.sum(v)
    return v

def forward_viterbi_wikipedia(obs, states, start_p, trans_p, emit_p):
    """
    This implementation is from wikipedia.
    It returns a triple.
    The first returned element is the probability of the observed sequence.
    The second returned element is the viterbi path.
    The third returned element is the probability of the viterbi path.
    Note that the returned viterbi path has a length that is one more than the observed sequence.
    To get a path that is conformant with the observed sequence,
    remove the last element of the viterbi path.
    @param obs: the observed sequence of emitted states
    @param states: list of all possible hidden states
    @param start_p: initial hidden state distribution as a dict
    @param trans_p: hidden state transition matrix as a dict of dicts
    @param emit_p: emission distributions as a dict of dicts
    @return: (total, argmax, valmax)
    """
    T = {}
    for state in states:
        T[state] = (start_p[state], [state], start_p[state])
    for output in obs:
        U = {}
        for next_state in states:
            total = 0
            argmax = None
            valmax = 0
            for source_state in states:
                (prob, v_path, v_prob) = T[source_state]
                p = emit_p[source_state][output] * trans_p[source_state][next_state]
                prob *= p
                v_prob *= p
                total += prob
                if v_prob > valmax:
                    argmax = v_path + [next_state]
                    valmax = v_prob
            U[next_state] = (total, argmax, valmax)
        T = U
    total = 0
    argmax = None
    valmax = 0
    for state in states:
        prob, v_path, v_prob = T[state]
        total += prob
        if v_prob > valmax:
            argmax = v_path
            valmax = v_prob
    return (total, argmax, valmax)


class TrainedModel:
    """
    This model includes transitions among hidden states.
    """
    def __init__(self, transition_matrix, hidden_state_objects):
        """
        The two arguments are conformantly ordered.
        @param transition_matrix: a numpy array that is a transition matrix among hidden states
        @param hidden_state_objects: a conformant list of hidden state objects
        """
        self.transition_matrix = transition_matrix
        self.hidden_state_objects = hidden_state_objects
        self.stationary_distribution = TransitionMatrix.get_stationary_distribution(transition_matrix)
        self.initial_distribution = self.stationary_distribution

    def set_initial_distribution(self, distribution):
        """
        The initial distribution might not be the stationary distribution.
        @param distribution: the desired initial distribution of hidden states
        """
        self.initial_distribution = distribution

    def sample_next_state_index(self, current_state_index):
        """
        @param current_state_index: the index of the current state
        @return: a state index sampled according to the transition matrix
        """
        return Util.random_weighted_int(self.transition_matrix[current_state_index])

    def get_joint_log_likelihood(self, hidden_seq, observed_seq):
        """
        The two arguments are conformantly ordered.
        @param hidden_seq: a sequence of hidden state indices
        @param observed_seq: a conformant sequence of observation objects
        @return: the joint likelihood of the hidden and observed sequences
        """
        # do validation
        if len(hidden_seq) != len(observed_seq):
            raise ValueError('expected conformant input sequences')
        # initialize the log likelihood
        log_accum = 0
        # add the contribution of the initial hidden state
        initial_hidden_state = hidden_seq[0]
        log_accum += math.log(self.initial_distribution[initial_hidden_state])
        # add the contribution of hidden state transitions
        for i, j in Util.pairwise(hidden_seq):
            log_accum += math.log(self.transition_matrix[i, j])
        # add the contribution of emissions
        for i, observation in zip(hidden_seq, observed_seq):
            log_accum += self.hidden_state_objects[i].get_log_likelihood(observation)
        # return the log likelihood
        return log_accum

    def forward_viterbi_wikipedia(self, observations):
        """
        @param observations: the sequence of observations
        @return: (observation probability, viterbi path, viterbi path probability)
        """
        nhidden = len(self.hidden_state_objects)
        T = {}
        for hidden_state in range(nhidden):
            p = self.initial_distribution[hidden_state]
            T[hidden_state] = (p, [hidden_state], p)
        for output in observations:
            U = {}
            for next_hidden_state in range(nhidden):
                total = 0
                argmax = None
                valmax = 0
                for source_hidden_state in range(nhidden):
                    prob, v_path, v_prob = T[source_hidden_state]
                    p = self.hidden_state_objects[source_hidden_state].get_likelihood(output)
                    p *= self.transition_matrix[source_hidden_state, next_hidden_state]
                    prob *= p
                    v_prob *= p
                    total += prob
                    if v_prob > valmax:
                        argmax = v_path + [next_hidden_state]
                        valmax = v_prob
                U[next_hidden_state] = (total, argmax, valmax)
            T = U
        total = 0
        argmax = None
        valmax = 0
        for hidden_state in range(nhidden):
            prob, v_path, v_prob = T[hidden_state]
            total += prob
            if v_prob > valmax:
                argmax = v_path
                valmax = v_prob
        return total, argmax, valmax

    def brute_posterior_decoding(self, observations):
        """
        Get the distribution of hidden states at each position given the observed sequence.
        This is done inefficiently by summing over each possible hidden state sequence.
        @param observations: the sequence of observations
        @return: hidden state distributions at each position, and total probability
        """
        nhidden = len(self.hidden_state_objects)
        total_log_sums = []
        # precalculate the log likelihood for each observation for each hidden state
        position_log_likelihoods = []
        for obs in observations:
            log_likelihoods = [state.get_log_likelihood(obs) for state in self.hidden_state_objects]
            position_log_likelihoods.append(log_likelihoods)
        # each hidden state at each position gets a list of log likelihoods
        total_accum = [[[] for i in range(nhidden)] for j in observations]
        # calculate the log likelihood for each hidden sequence
        for hidden_sequence in itertools.product(range(nhidden), repeat=len(observations)):
            accum = 0
            accum += math.log(self.initial_distribution[hidden_sequence[0]])
            for i, j in Util.pairwise(hidden_sequence):
                accum += math.log(self.transition_matrix[i, j])
            for index, log_likelihoods in zip(hidden_sequence, position_log_likelihoods):
                accum += log_likelihoods[index]
            # accumulate the log likelihood
            for i, hidden_state in enumerate(hidden_sequence):
                total_accum[i][hidden_state].append(accum)
            # add to the total probability
            total_log_sums.append(accum)
        # get the distribution at each position
        distributions = []
        for log_distribution_lists in total_accum:
            distribution = [scipy.maxentropy.logsumexp(x) for x in log_distribution_lists]
            distribution = [d - max(distribution) for d in distribution]
            distribution = [math.exp(d) for d in distribution]
            distribution = [d / sum(distribution) for d in distribution]
            distributions.append(distribution)
        total_probability = math.exp(scipy.maxentropy.logsumexp(total_log_sums))
        return distributions, total_probability

    def naive_forward_durbin(self, observations):
        """
        Implement the forward algorithm directly from the book.
        The book is Biological Sequence Analysis by Durbin et al.
        @param observations: the sequence of observations
        @return: the list over positions of the list of forward variables over states, and total probability
        """
        nhidden = len(self.hidden_state_objects)
        nobs = len(observations)
        f = [[0]*nhidden for i in range(nobs)]
        # define the initial f variable
        for sink_index, sink_state in enumerate(self.hidden_state_objects):
            f[0][sink_index] = sink_state.get_likelihood(observations[0]) * self.initial_distribution[sink_index]
        # define the subsequent f variables
        for i in range(1, nobs):
            obs = observations[i]
            for sink_index, sink_state in enumerate(self.hidden_state_objects):
                f[i][sink_index] = sink_state.get_likelihood(obs)
                p = 0
                for source_index, source_state in enumerate(self.hidden_state_objects):
                    p += f[i-1][source_index] * self.transition_matrix[source_index, sink_index]
                f[i][sink_index] *= p
        total_probability = 0
        for source_index, source_state in enumerate(self.hidden_state_objects):
            total_probability += f[nobs-1][source_index]
        return f, total_probability

    def naive_backward_durbin(self, observations):
        """
        Implement the backward algorithm directly from the book.
        The book is Biological Sequence Analysis by Durbin et al.
        @param observations: the sequence of observations
        @return: the list over positions of the list of backward variables over states, and total probability
        """
        nhidden = len(self.hidden_state_objects)
        nobs = len(observations)
        b = [[0]*nhidden for i in observations]
        b[nobs-1] = [1]*nhidden
        for i in reversed(range(nobs-1)):
            for source_index, source_state in enumerate(self.hidden_state_objects):
                for sink_index, sink_state in enumerate(self.hidden_state_objects):
                    p = 1.0
                    p *= self.transition_matrix[source_index, sink_index]
                    p *= sink_state.get_likelihood(observations[i+1])
                    p *= b[i+1][sink_index]
                    b[i][source_index] += p
        total_probability = 0
        for sink_index, sink_state in enumerate(self.hidden_state_objects):
            p = 1.0
            p *= self.initial_distribution[sink_index]
            p *= sink_state.get_likelihood(observations[0])
            p *= b[0][sink_index]
            total_probability += p
        return b, total_probability

    def naive_posterior_durbin(self, observations):
        """
        Implement the backward algorithm directly from the book.
        The book is Biological Sequence Analysis by Durbin et al.
        @param observations: the sequence of observations
        @return: the list over positions of the posterior hidden state distributions, and the total probability
        """
        f, total_f = self.naive_forward_durbin(observations)
        b, total_b = self.naive_backward_durbin(observations)
        if not np.allclose(total_f, total_b):
            raise ValueError('inconsistent total probability calculations')
        total = (total_f + total_b) / 2
        distributions = []
        for fs, bs in zip(f, b):
            distribution = [x*y/total for x, y in zip(fs, bs)]
            if not np.allclose(sum(distribution), 1):
                raise ValueError('the distribution does not sum to 1: ' + str(sum(distribution)))
            distributions.append(distribution)
        return distributions, total

    def scaled_forward_durbin(self, observations):
        """
        Implement the scaled forward algorithm directly from the book.
        The book is Biological Sequence Analysis by Durbin et al.
        At each position, the sum over states of the f variable is 1.
        @param observations: the sequence of observations
        @return: the list of lists of scaled f variables, and the scaling variables
        """
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
            # define an unscaled f variable at this position
            for sink_index, sink_state in enumerate(self.hidden_state_objects):
                f[i][sink_index] = sink_state.get_likelihood(obs)
                p = 0
                for source_index, source_state in enumerate(self.hidden_state_objects):
                    p += f[i-1][source_index] * self.transition_matrix[source_index, sink_index]
                f[i][sink_index] *= p
            # define the scaling factor at this position
            s[i] = sum(f[i])
            # define the scaled f variable at this position
            for sink_index in range(nhidden):
                f[i][sink_index] /= s[i]
        return f, s

    def scaled_backward_durbin(self, observations, scaling_factors):
        """
        Implement the scaled backward algorithm directly from the book.
        The book is Biological Sequence Analysis by Durbin et al.
        The scaling factors must have been calculated using the scaled forward algorithm.
        @param observations: the sequence of observations
        @param scaling_factors: the scaling factor for each position
        @return: the list of lists of scaled b variables
        """
        nhidden = len(self.hidden_state_objects)
        nobs = len(observations)
        b = [[0]*nhidden for i in observations]
        b[nobs-1] = [1/scaling_factors[nobs-1]]*nhidden
        for i in reversed(range(nobs-1)):
            for source_index in range(nhidden):
                accum = 0
                for sink_index, sink_state in enumerate(self.hidden_state_objects):
                    p = 1.0
                    p *= self.transition_matrix[source_index, sink_index]
                    p *= sink_state.get_likelihood(observations[i+1])
                    p *= b[i+1][sink_index]
                    accum += p
                b[i][source_index] = accum / scaling_factors[i]
        return b

    def scaled_posterior_durbin(self, observations):
        """
        Implement the backward algorithm directly from the book.
        The book is Biological Sequence Analysis by Durbin et al.
        @param observations: the sequence of observations
        @return: the list of position-specific posterior hidden state distributions
        """
        f, s = self.scaled_forward_durbin(observations)
        b = self.scaled_backward_durbin(observations, s)
        distributions = []
        #for i, (fs, bs) in enumerate(zip(f, b)):
        for fs, bs, si in zip(f, b, s):
            distribution = [x*y*si for x, y in zip(fs, bs)]
            distributions.append(distribution)
        return distributions

    def naive_transition_expectations_durbin(self, observations):
        """
        This is part of the Baum-Welch transition matrix estimation.
        To get the estimated transition matrix,
        normalize the matrix of expected transitions.
        @param observations: the sequence of observations
        @return: expected initial counts, and expected transition counts
        """
        nhidden = len(self.hidden_state_objects)
        nobs = len(observations)
        # get the naive (unscaled) forward and backward variables and the total probability
        f, total_f = self.naive_forward_durbin(observations)
        b, total_b = self.naive_backward_durbin(observations)
        if not np.allclose(total_f, total_b):
            raise ValueError('inconsistent total probability calculations')
        total = (total_f + total_b) / 2
        # initialize the matrix of expected initial states
        initial_counts = [0]*nhidden
        # initialize the matrix of expected counts
        A = np.zeros((nhidden, nhidden))
        # get the expected counts for each transition
        for sink_index, sink_state in enumerate(self.hidden_state_objects):
            # get expected initial state counts
            p = 1.0
            p *= self.initial_distribution[sink_index]
            p *= sink_state.get_likelihood(observations[0])
            p *= b[0][sink_index]
            initial_counts[sink_index] = p / total
            # get expected transition counts
            for source_index, source_state in enumerate(self.hidden_state_objects):
                for i in range(1, nobs):
                    obs = observations[i]
                    p = 1.0
                    p *= f[i-1][source_index]
                    p *= self.transition_matrix[source_index, sink_index]
                    p *= sink_state.get_likelihood(obs)
                    p *= b[i][sink_index]
                    A[source_index, sink_index] += p
                A[source_index, sink_index] /= total
        return initial_counts, A

    def scaled_transition_expectations_durbin(self, observations):
        """
        This is part of the Baum-Welch transition matrix estimation.
        To get the estimated transition matrix,
        normalize the matrix of expected transitions.
        @param observations: the sequence of observations
        @return: expected initial counts, and expected transition counts
        """
        nhidden = len(self.hidden_state_objects)
        nobs = len(observations)
        # get the scaled forward and backward variables
        f, s = self.scaled_forward_durbin(observations)
        b = self.scaled_backward_durbin(observations, s)
        # initialize the matrix of expected initial states
        initial_counts = [0]*nhidden
        # initialize the matrix of expected counts
        A = np.zeros((nhidden, nhidden))
        # get the expected counts for each transition
        for sink_index, sink_state in enumerate(self.hidden_state_objects):
            # get expected initial state counts
            p = 1.0
            p *= self.initial_distribution[sink_index]
            p *= sink_state.get_likelihood(observations[0])
            p *= b[0][sink_index]
            initial_counts[sink_index] = p
            # get expected transition counts
            for source_index, source_state in enumerate(self.hidden_state_objects):
                for i in range(1, nobs):
                    obs = observations[i]
                    p = 1.0
                    p *= f[i-1][source_index]
                    p *= self.transition_matrix[source_index, sink_index]
                    p *= sink_state.get_likelihood(obs)
                    p *= b[i][sink_index]
                    A[source_index, sink_index] += p
        return initial_counts, A

    def naive_emission_expectations_durbin(self, observations):
        """
        Get the expected emission distribution for each state.
        With small modifications this function can be used to do more useful things.
        @param observations: the sequence of observations
        @return: for each hidden state, the expected count of each emitted state
        """
        nhidden = len(self.hidden_state_objects)
        nobs = len(observations)
        # get the naive (unscaled) forward and backward variables and the total probability
        f, total_f = self.naive_forward_durbin(observations)
        b, total_b = self.naive_backward_durbin(observations)
        if not np.allclose(total_f, total_b):
            raise ValueError('inconsistent total probability calculations')
        total = (total_f + total_b) / 2
        # initialize the expected emission counts for each observed state
        E = [{} for i in range(nhidden)]
        # count the expected number of each emitted observation state
        for i, obs in enumerate(observations):
            for source_index in range(nhidden):
                count = f[i][source_index] * b[i][source_index] / total
                E[source_index][obs] = E[source_index].get(obs, 0) + count
        return E

    def scaled_emission_expectations_durbin(self, observations):
        """
        Get the expected emission distribution for each state.
        With small modifications this function can be used to do more useful things.
        @param observations: the sequence of observations
        @return: for each hidden state, the expected count of each emitted state
        """
        nhidden = len(self.hidden_state_objects)
        nobs = len(observations)
        # get the scaled forward and backward variables and the scaling factors
        f, s = self.scaled_forward_durbin(observations)
        b = self.scaled_backward_durbin(observations, s)
        # initialize the expected emission counts for each observed state
        E = [{} for i in range(nhidden)]
        # count the expected number of each emitted observation state
        for i, obs in enumerate(observations):
            for source_index in range(nhidden):
                count = f[i][source_index] * b[i][source_index] * s[i]
                E[source_index][obs] = E[source_index].get(obs, 0) + count
        return E

class DishonestCasino(TrainedModel):
    def __init__(self):
        fair_state = HiddenDieState(1/6.0)
        loaded_state = HiddenDieState(0.5)
        transition_matrix = np.array([[0.95, 0.05], [0.1, 0.9]])
        TrainedModel.__init__(self, transition_matrix, [fair_state, loaded_state])

class CrazyCasino(TrainedModel):
    def __init__(self):
        fair_state = HiddenDieState(1/6.0)
        loaded_state = OnlySixes()
        transition_matrix = np.array([[0.95, 0.05], [0.1, 0.9]])
        TrainedModel.__init__(self, transition_matrix, [fair_state, loaded_state])


class HiddenState:
    """
    Each hidden state emits some object according to a distribution.
    The hidden state object has two main roles.
    Firstly, it allows an observation to be sampled.
    Secondly, it computes the likelihood of an observation.
    This is a virtual base class.
    """
    def sample_observation(self):
        raise NotImplementedError()
    def get_likelihood(self, observation):
        raise NotImplementedError()
    def get_log_likelihood(self, observation):
        raise NotImplementedError()


class HiddenDieState:
    """
    This is for the dishonest casino example.
    """
    def __init__(self, probability_of_six):
        if not (0 <= probability_of_six <= 1):
            raise ValueError('probability_of_six must be between 0 and 1')
        self.probability_of_six = probability_of_six

    def sample_observation(self):
        if random.random() < self.probability_of_six:
            return 6
        else:
            return random.choice((1, 2, 3, 4, 5))

    def get_likelihood(self, observation):
        if observation not in ((1, 2, 3, 4, 5, 6)):
            raise ValueError()
        elif observation == 6:
            return self.probability_of_six
        else:
            return (1.0 - self.probability_of_six) / 5.0

    def get_log_likelihood(self, observation):
        return math.log(self.get_likelihood(observation))


class OnlySixes:
    def sample_observation(self):
        return 6
    def get_likelihood(self, observation):
        return 1.0 if observation == 6 else 0.0
    def get_log_likelihood(self, observation):
        return 0.0 if observation == 6 else float('-inf')


class TestHMM(unittest.TestCase):

    def assert_almost_equal_emissions(self, Ea, Eb):
        """
        This is a helper function.
        Assert that these lists of emission dictionaries are almost the same.
        @param Ea: a list of emission dictionaries
        @param Eb: a list of emission dictionaries
        """
        for distribution_a, distribution_b in zip(Ea, Eb):
            # assert that each dictionary has the same set of keys
            keys_a = set(distribution_a.keys())
            keys_b = set(distribution_b.keys())
            self.assertEqual(keys_a, keys_b)
            # assert that for each key the entries are almost the same
            for k in keys_a:
                self.assertAlmostEqual(distribution_a[k], distribution_b[k])

    def test_dishonest_casino_sampling(self):
        state = HiddenDieState(.5)
        observations = [state.sample_observation() for i in range(10)]
        likelihoods = [state.get_likelihood(x) for x in observations]

    def test_dishonest_casino_viterbi_wikipedia(self):
        """
        Test the occasionally dishonest casino HMM with the wikipedia viterbi code.
        Use the code with nested dictionaries.
        """
        hidden_states = ['F', 'L']
        start_p = {'F' : 0.5, 'L' : 0.5}
        trans_p = {
                'F' : {'F' : 0.95, 'L' : 0.05},
                'L' : {'F' : 0.1, 'L' : 0.9}}
        emit_p = {
                'F' : {1 : 1/6., 2 : 1/6., 3 : 1/6., 4 : 1/6., 5 : 1/6., 6 : 1/6.},
                'L' : {1 : 0.1, 2 : 0.1, 3 : 0.1, 4 : 0.1, 5 : 0.1, 6 : 0.5}}
        observations, estimates = get_example_rolls()
        expected_viterbi_path = ''.join(estimates)
        total, argmax, valmax = forward_viterbi_wikipedia(observations, hidden_states, start_p, trans_p, emit_p)
        observed_viterbi_path = ''.join(argmax)[:-1]
        self.assertEqual(expected_viterbi_path, observed_viterbi_path)

    def test_dishonest_casino_viterbi_objects(self):
        """
        Test the occasionally dishonest casino HMM with the wikipedia viterbi code.
        Use the code with hidden markov model objects.
        """
        hmm = DishonestCasino()
        # define the expected result
        observations, estimates = get_example_rolls()
        expected_viterbi_path = ''.join(estimates)
        # define the result we actually get
        total, argmax, valmax = hmm.forward_viterbi_wikipedia(observations)
        state_to_name = ['F', 'L']
        observed_viterbi_path = ''.join(state_to_name[i] for i in argmax[:-1])
        # assert that the results are the same
        self.assertEqual(expected_viterbi_path, observed_viterbi_path)

    def test_brute_posterior_decoding(self):
        hmm = DishonestCasino()
        # define a sequence of observations
        observations = [1, 2, 6, 6, 1, 2, 3, 4, 5, 6]
        # get the posterior hidden marginal distribution at each position
        wikipedia_total, argmax, valmax = hmm.forward_viterbi_wikipedia(observations)
        distributions, brute_total = hmm.brute_posterior_decoding(observations)
        # assert that the totals are the same
        self.assertAlmostEqual(wikipedia_total, brute_total)
        # assert that the distributions are reasonable
        for distribution in distributions:
            self.assertAlmostEqual(sum(distribution), 1)

    def test_naive_posterior_decoding(self):
        hmm = DishonestCasino()
        # define a sequence of observations
        observations = [1, 2, 6, 6, 1, 2, 3, 4, 5, 6]
        # get the posterior hidden marginal distribution at each position
        brute_distributions, total_brute = hmm.brute_posterior_decoding(observations)
        naive_distributions, total_naive = hmm.naive_posterior_durbin(observations)
        f, total_f = hmm.naive_forward_durbin(observations)
        b, total_b = hmm.naive_backward_durbin(observations)
        # assert that the totals are the same
        self.assertAlmostEqual(total_brute, total_naive)
        self.assertAlmostEqual(total_brute, total_f)
        self.assertAlmostEqual(total_brute, total_b)
        # assert that the posterior distributions are the same
        for bdist, ndist in zip(brute_distributions, naive_distributions):
            self.assertTrue(np.allclose(bdist, ndist))

    def test_scaled_forward_durbin(self):
        hmm = DishonestCasino()
        # define a sequence of observations
        observations = [1, 2, 6, 6, 1, 2, 3, 4, 5, 6]
        # Get the forward variables and total probability
        # using the unscaled forward algorithm.
        f, total_probability = hmm.naive_forward_durbin(observations)
        # Get the scaled forward variables and scaling factors
        # using the scaled forward algorithm.
        f, s = hmm.scaled_forward_durbin(observations)
        # assert that the product of the scaling factors is the total probability
        self.assertAlmostEqual(Util.product(s), total_probability)

    def test_scaled_posterior_durbin(self):
        hmm = DishonestCasino()
        # define a sequence of observations
        observations = [1, 2, 6, 6, 1, 2, 3, 4, 5, 6]
        # Get posterior distributions using a naive algorithm.
        naive_distributions, total_naive = hmm.naive_posterior_durbin(observations)
        # Get posterior distributions using a scaled algorithm.
        scaled_distributions = hmm.scaled_posterior_durbin(observations)
        # assert that the distributions are the same
        for dnaive, dscaled in zip(naive_distributions, scaled_distributions):
            self.assertTrue(np.allclose(dnaive, dscaled))

    def test_naive_transition_expectations_durbin(self):
        hmm = DishonestCasino()
        # define a sequence of observations
        observations = [1, 2, 6, 6, 1, 2, 3, 4, 5, 6]
        # Get expected initial state counts and expected transition counts.
        initial_counts, A = hmm.naive_transition_expectations_durbin(observations)
        # the initial counts should sum to 1 because there is one initial state
        self.assertAlmostEqual(sum(initial_counts), 1)
        # the sum of elements of A should be one fewer than the number of observations
        self.assertAlmostEqual(np.sum(A), len(observations)-1)

    def test_scaled_transition_expectations_durbin(self):
        hmm = DishonestCasino()
        # define a sequence of observations
        observations = [1, 2, 6, 6, 1, 2, 3, 4, 5, 6]
        # Get expected initial state counts and expected transition counts.
        naive_initial, naive_A = hmm.naive_transition_expectations_durbin(observations)
        scaled_initial, scaled_A = hmm.scaled_transition_expectations_durbin(observations)
        # the results should be the same for each method
        self.assertTrue(np.allclose(naive_initial, scaled_initial))
        self.assertTrue(np.allclose(naive_A, scaled_A))

    def test_naive_emission_expectations_durbin(self):
        hmm = DishonestCasino()
        # define a sequence of observations
        observations = [1, 2, 6, 6, 1, 2, 3, 4, 5, 6]
        # get emission count expectations
        E = hmm.naive_emission_expectations_durbin(observations)
        # get the emission count expectation sums per hidden state
        expectation_per_hidden_state = [sum(d.values()) for d in E]
        # Get the posterior distribution per position,
        # and then sum over positions to get the posterior
        # distribution of hidden states.
        hidden_distribution_per_position = hmm.scaled_posterior_durbin(observations)
        hidden_counts = [sum(d) for d in zip(*hidden_distribution_per_position)]
        # assert that the two count vectors are the same
        self.assertTrue(np.allclose(expectation_per_hidden_state, hidden_counts))

    def test_scaled_emission_expectations_durbin(self):
        hmm = DishonestCasino()
        # define a sequence of observations
        observations = [1, 2, 6, 6, 1, 2, 3, 4, 5, 6]
        # get emission count expectations using naive and scaled algorithms
        E_scaled = hmm.scaled_emission_expectations_durbin(observations)
        E_naive = hmm.naive_emission_expectations_durbin(observations)
        # assert that they are almost the same
        self.assert_almost_equal_emissions(E_scaled, E_naive)

    def test_emission_posterior_equivalence(self):
        """
        See if two things are the same.
        If you do posterior decoding,
        and then get the distribution of observations per hidden state,
        see if this is the same as getting the emission distribution directly.
        """
        hmm = DishonestCasino()
        nhidden = 2
        # define a sequence of observations
        observations = [1, 2, 6, 6, 1, 2, 3, 4, 5, 6]
        # get emission count expectations using the scaled algorithm
        E_scaled = hmm.scaled_emission_expectations_durbin(observations)
        # get the posterior distribution per position
        posterior_per_position = hmm.scaled_posterior_durbin(observations)
        # for each hidden state get a posterior expectation of counts
        E_hacked = [{} for i in range(nhidden)]
        for obs, posterior in zip(observations, posterior_per_position):
            for hidden_state in range(nhidden):
                count = posterior[hidden_state]
                E_hacked[hidden_state][obs] = E_hacked[hidden_state].get(obs, 0) + count
        # assert that they are almost the same
        self.assert_almost_equal_emissions(E_scaled, E_hacked)

    def test_crazy_casino(self):
        """
        Test robustness vs. degenerate emission distributions.
        One of the hidden states has zero probability of emitting states
        that are present in the observed sequence.
        """
        hmm = CrazyCasino()
        observations = [1, 2, 6, 6, 1, 2, 3, 4, 5, 6]
        posterior_per_position = hmm.scaled_posterior_durbin(observations)
        # each posterior distribution should sum to 1
        for posterior in posterior_per_position:
            self.assertAlmostEqual(sum(posterior), 1)
        # the probability of fairness at the second position should be 1
        self.assertAlmostEqual(posterior_per_position[1][0], 1)
        # neither probability at the third position should be almost zero
        self.assertNotAlmostEqual(posterior_per_position[2][0], 0)
        self.assertNotAlmostEqual(posterior_per_position[2][1], 0)

    def test_uninformative_transition_matrix(self):
        """
        Test an uninformative transition matrix.
        When the transition matrix is uninformative,
        the position specific hidden state posterior distribution
        can be calculated without using a Markov model.
        """
        # define the model
        nhidden = 2
        T = np.ones((nhidden, nhidden)) / float(nhidden)
        fair_state = HiddenDieState(1/6.0)
        loaded_state = HiddenDieState(0.5)
        models = [fair_state, loaded_state]
        hmm = TrainedModel(T, models)
        # define the observations
        observations = [1, 2, 6, 6, 1, 2, 3, 4, 5, 6]
        # get the posterior distributions using the hmm
        hmm_distributions = hmm.scaled_posterior_durbin(observations)
        # get the posterior distributions using only the individual models
        individual_distributions = [get_distribution(models, obs) for obs in observations]
        # assert that the distributions are the same
        self.assertTrue(np.allclose(hmm_distributions, individual_distributions))


if __name__ == '__main__':
    unittest.main()
