"""
Sample events along a path.
Do the sampling conditional on the beginning and ending states.
Assume a Markov process with a known rate matrix.
"""

import unittest
import math
import random

import numpy as np
import scipy

import Util
import RateMatrix


class MatrixPowerCache:
    """
    Get various powers of a matrix with caching for speed.
    """
    def __init__(self, matrix):
        self.matrix_powers = []
        self.matrix_powers.append(np.eye(matrix.shape[0], matrix.shape[1]))
        self.matrix_powers.append(matrix)

    def get_power(self, power):
        while len(self.matrix_powers) <= power:
            matrix = self.matrix_powers[1]
            last_matrix = self.matrix_powers[-1]
            new_matrix = np.dot(matrix, last_matrix)
            self.matrix_powers.append(new_matrix)
        return self.matrix_powers[power]



def get_discrete_path_sample(initial_state, terminal_state, states, path_length, transition_matrix):
    """
    Return a sequence of states starting at the initial state and ending at the terminal state.
    Consecutive states may be the same if the transition matrix has positive diagonal elements.
    @param initial_state: the first state
    @param terminal_state: the last state
    @param states: the ordered names of the states
    @param path_length: the number of states in the returned path
    @type path_length: int
    @param transition_matrix: a dictionary mapping an ordered state pair to a probability
    @type transition_matrix: dict
    """
    # assert that the states are valid
    state_set = set(states)
    assert len(state_set) == len(states)
    assert initial_state in state_set
    assert terminal_state in state_set
    for a, b in transition_matrix:
        assert a in state_set
        assert b in state_set
    assert path_length == int(path_length)
    # take care of degenerate cases
    if path_length == 0:
        return []
    elif path_length == 1:
        assert initial_state == terminal_state
        return [initial_state]
    elif path_length == 2:
        return [initial_state, terminal_state]
    # create transition matrices raised to various powers
    max_power = path_length - 2
    matrix_powers = [1]
    state_to_index = dict((state, i) for i, state in enumerate(states))
    m = np.zeros((len(states), len(states)))
    for (a, b), p in transition_matrix.items():
        ia = state_to_index[a]
        ib = state_to_index[b]
        m[ia, ib] = p
    matrix_powers.append(m)
    for i in range(2, max_power + 1):
        matrix_powers.append(np.dot(matrix_powers[i-1], m))
    # sample the path
    path = [initial_state]
    for i in range(1, path_length-1):
        previous_state = path[i-1]
        weight_state_pairs = []
        for state_index, state in enumerate(states):
            weight = 1
            weight *= m[state_to_index[previous_state], state_index]
            weight *= matrix_powers[path_length - 1 - i][state_index, state_to_index[terminal_state]]
            weight_state_pairs.append((weight, state))
        next_state = Util.weighted_choice(weight_state_pairs)
        path.append(next_state)
    path.append(terminal_state)
    return path


def get_uniformization_sample(initial_state, terminal_state, states, path_length, rate_matrix):
    """
    Follow the explanation in a manuscript by Stone and Hobolth.
    scipy is a non-standard python library
    @return: a list of (time, state) state change events
    """
    # map states to indices
    state_to_index = dict((state, i) for i, state in enumerate(states))
    # find the maximum rate away from a state
    max_rate = max(-rate_matrix[(a, a)] for a in states)
    # create a uniformized discrete transition matrix in convenient dictionary form
    discrete_transition_matrix = {}
    for (a, b), r in rate_matrix.items():
        discrete_transition_matrix[(a, b)] = r / max_rate
        if a == b:
            discrete_transition_matrix[(a, b)] += 1.0
    # create a discrete transition matrix in the numpy format,
    # and create the rate matrix in the numpy format
    R = np.zeros((len(states), len(states)))
    numpy_rate_matrix = np.zeros((len(states), len(states)))
    for (a, b), r in rate_matrix.items():
        ia = state_to_index[a]
        ib = state_to_index[b]
        numpy_rate_matrix[ia, ib] = r
        R[ia, ib] = discrete_transition_matrix[(a, b)]
    # convert initial and terminal states to indices
    initial_index = state_to_index[initial_state]
    terminal_index = state_to_index[terminal_state]
    # get the probability of the terminal state given the initial state and the path length
    rate_matrix_exponential = scipy.linalg.matfuncs.expm(numpy_rate_matrix * path_length)
    Pab = rate_matrix_exponential[initial_index, terminal_index]
    # draw the number of state changes
    cumulative_probability = 0
    n = 0
    matrix_powers = MatrixPowerCache(R)
    cutoff = random.uniform(0, Pab)
    #print 'cutoff =', cutoff
    #print 'initial_index =', initial_index
    #print 'terminal_index =', terminal_index
    #print matrix_powers.get_power(0)
    while 1:
        poisson_factor = scipy.stats.poisson.pmf(n, max_rate * path_length)
        discrete_transition_factor = matrix_powers.get_power(n)[initial_index, terminal_index]
        cumulative_probability += poisson_factor * discrete_transition_factor
        #print 'cumulative probability =', cumulative_probability
        if cutoff < cumulative_probability:
            break
        n += 1
    #print 'n =', n
    # deal with degenerate cases
    if n == 0:
        return []
    elif n == 1:
        if initial_state == terminal_state:
            return []
        else:
            return [(random.uniform(0, path_length), terminal_state)]
    # Simulate a discrete path given the number of changes and the initial and terminal states.
    # The path is called virtual because some changes may be from a state to itself.
    virtual_path = get_discrete_path_sample(initial_state, terminal_state, states, n+1, discrete_transition_matrix)[1:]
    virtual_times = list(sorted(random.uniform(0, path_length) for i in range(n)))
    events = []
    last_state = initial_state
    last_time = 0
    for current_state, current_time in zip(virtual_path, virtual_times):
        if current_state == last_state:
            continue
        events.append((current_state, current_time))
        last_state = current_state
        last_time = current_time
    return events


def get_nielsen_sample(initial_state, terminal_state, states, path_length, rate_matrix):
    """
    Inspired by a 2001 paper by Rasmus Nielsen.
    Genetics. 2001 September; 159(1): 401-411.
    Mutations as missing data: inferences on the ages and distributions of nonsynonymous and synonymous mutations.
    @param initial_state: the initial state of the path
    @param terminal_state: the terminal state of the path
    @param states: the states allowed in the model
    @param path_length: the length or time of the path, depending on its interpretation
    @param rate_matrix: this object defines the Markov process of the path
    @return: a list of (time, state) events or None if the attempt was rejected
    """
    # sample events
    t = 0
    events = []
    state = initial_state
    if terminal_state != initial_state:
        rate = -rate_matrix[(state, state)]
        U = random.uniform(0, 1)
        t = -math.log(1.0 - U*(1.0 - math.exp(-path_length*rate)))/rate
        assert t < path_length
        weight_state_pairs = [(rate_matrix[(state, next)], next) for next in states if next != state]
        state = Util.weighted_choice(weight_state_pairs)
        events.append((t, state))
    while True:
        rate = -rate_matrix[(state, state)]
        t += random.expovariate(rate)
        if t > path_length:
            if state == terminal_state:
                return events
            else:
                return None
        weight_state_pairs = [(rate_matrix[(state, next)], next) for next in states if next != state]
        state = Util.weighted_choice(weight_state_pairs)
        events.append((t, state))
    return events

def get_rejection_sample(initial_state, terminal_state, states, path_length, rate_matrix):
    """
    @param initial_state: the initial state of the path
    @param terminal_state: the terminal state of the path
    @param states: the states allowed in the model
    @param path_length: the length or time of the path, depending on its interpretation
    @param rate_matrix: this object defines the Markov process of the path
    @return: a list of (time, state) events or None if the attempt was rejected
    """
    events = []
    t = 0
    state = initial_state
    while True:
        rate_away = -rate_matrix[(state, state)]
        t += random.expovariate(rate_away)
        if t > path_length:
            if state == terminal_state:
                return events
            else:
                return None
        weight_state_pairs = [(rate_matrix[(state, next)], next) for next in states if next != state]
        state = Util.weighted_choice(weight_state_pairs)
        events.append((t, state))


class TestPathSampler(unittest.TestCase):

    def test_jukes_cantor_rejection(self):
        path_length = 1
        jukes_cantor_rate_matrix = RateMatrix.get_jukes_cantor_rate_matrix()
        states = 'ACGT'
        n = 200
        observed = 0
        for i in range(n):
            events = get_rejection_sample('A', 'C', states, path_length, jukes_cantor_rate_matrix)
            if events is not None:
                observed += 1
        p = RateMatrix.get_jukes_cantor_transition_matrix(path_length)[('A', 'C')]
        expected = n*p
        variance = n*p*(1-p)
        errstr = 'observed: %f  expected: %f' % (observed, expected)
        self.failUnless(abs(observed - expected) < 3*math.sqrt(variance), errstr)

    def test_hky_nielsen(self):
        """
        Give modified rejection sampling a chance to fail.
        It should give the same results as vanilla rejection sampling.
        """
        distribution = {'A':.2,'C':.3,'G':.3,'T':.2}
        kappa = 2
        rate_matrix_object = RateMatrix.get_unscaled_hky85_rate_matrix(distribution, kappa)
        rate_matrix_object.normalize()
        rate_matrix = rate_matrix_object.get_dictionary_rate_matrix()
        path_length = 2
        initial_state = 'A'
        terminal_state = 'C'
        states = 'ACGT'
        iterations = 200
        rejection_changes = []
        i = 0
        while i < iterations:
            rejection_events = get_rejection_sample(initial_state, terminal_state, states, path_length, rate_matrix)
            if rejection_events is not None:
                rejection_changes.append(len(rejection_events))
                i += 1
        nielsen_changes = []
        i = 0
        while i < iterations:
            nielsen_events = get_nielsen_sample(initial_state, terminal_state, states, path_length, rate_matrix)
            if nielsen_events is not None:
                nielsen_changes.append(len(nielsen_events))
                i += 1
        t, p = scipy.stats.mannwhitneyu(rejection_changes, nielsen_changes)
        self.failIf(p < .001)

    def test_hky_uniformization(self):
        """
        Give uniformization a chance to fail.
        It should give the same results as modified rejection sampling.
        """
        distribution = {'A':.2,'C':.3,'G':.3,'T':.2}
        kappa = 2
        rate_matrix_object = RateMatrix.get_unscaled_hky85_rate_matrix(distribution, kappa)
        rate_matrix_object.normalize()
        rate_matrix = rate_matrix_object.get_dictionary_rate_matrix()
        path_length = 2
        initial_state = 'A'
        terminal_state = 'C'
        states = 'ACGT'
        iterations = 200
        # get the modified rejection sampling changes, where each change is the number of events on a sampled path
        nielsen_changes = []
        i = 0
        while i < iterations:
            nielsen_events = get_nielsen_sample(initial_state, terminal_state, states, path_length, rate_matrix)
            if nielsen_events is not None:
                nielsen_changes.append(len(nielsen_events))
                i += 1
        # get the uniformization changes, where each change is the number of events on a sampled path
        uniformization_changes = []
        for i in range(iterations):
            uniformization_events = get_uniformization_sample(initial_state, terminal_state, states, path_length, rate_matrix)
            uniformization_changes.append(len(uniformization_events))
        # see if there is a statistically significant difference between the sampled path lengths
        #print sum(nielsen_changes)
        #print sum(uniformization_changes)
        t, p = scipy.stats.mannwhitneyu(uniformization_changes, nielsen_changes)
        self.failIf(p < .001, p)


def demo_rejection_sampling():
    path_length = 2
    jukes_cantor_rate_matrix = RateMatrix.get_jukes_cantor_rate_matrix()
    states = 'ACGT'
    n = 100000
    nielsen_event_count = 0
    nielsen_path_count = 0
    nielsen_first_time_sum = 0
    nielsen_dwell = dict((c, 0) for c in states)
    rejection_event_count = 0
    rejection_path_count = 0
    rejection_first_time_sum = 0
    rejection_dwell = dict((c, 0) for c in states)
    for i in range(n):
        initial_state = 'A'
        terminal_state = 'C'
        events = get_rejection_sample(initial_state, terminal_state, states, path_length, jukes_cantor_rate_matrix)
        if events is not None:
            assert events
            rejection_path_count += 1
            rejection_event_count += len(events)
            t, state = events[0]
            rejection_first_time_sum += t
            extended = [(0, initial_state)] + events + [(path_length, terminal_state)]
            for (t0, state0), (t1, state1) in zip(extended[:-1], extended[1:]):
                rejection_dwell[state0] += t1 - t0
        events = get_nielsen_sample(initial_state, terminal_state, states, path_length, jukes_cantor_rate_matrix)
        if events is not None:
            assert events
            nielsen_path_count += 1
            nielsen_event_count += len(events)
            t, state = events[0]
            nielsen_first_time_sum += t
            extended = [(0, initial_state)] + events + [(path_length, terminal_state)]
            for (t0, state0), (t1, state1) in zip(extended[:-1], extended[1:]):
                nielsen_dwell[state0] += t1 - t0
    expected_fraction = RateMatrix.get_jukes_cantor_transition_matrix(path_length)[(initial_state, terminal_state)]
    print 'testing the rejection sampling:'
    print 'expected fraction:', expected_fraction
    print 'observed fraction:', rejection_path_count / float(n)
    print 'comparing rejection sampling and nielsen sampling:'
    rejection_method_fraction = rejection_event_count / float(rejection_path_count)
    nielsen_method_fraction = nielsen_event_count / float(nielsen_path_count)
    print 'rejection method fraction:', rejection_method_fraction
    print 'nielsen method fraction:', nielsen_method_fraction
    print 'comparing time of first event:'
    print 'rejection method first event time mean:', rejection_first_time_sum / float(rejection_path_count)
    print 'nielsen method first event time mean:', nielsen_first_time_sum / float(nielsen_path_count)
    print 'comparing the duration spent in each state:'
    print 'rejection:'
    for state, t in rejection_dwell.items():
        print '\t%s: %f' % (state, t/float(rejection_path_count))
    print 'nielsen:'
    for state, t in nielsen_dwell.items():
        print '\t%s: %f' % (state, t/float(nielsen_path_count))

def demo_discrete_path_sampling():
    import EnglishModel
    initial_state = 'a'
    terminal_state = 'b'
    states = 'abcdefghijklmnopqrstuvwxyz '
    path_length = 10
    transition_matrix = EnglishModel.get_transition_matrix()
    path = get_discrete_path_sample(initial_state, terminal_state, states, path_length, transition_matrix)
    print ''.join(path)

def demo_uniformization():
    distribution = {'A':.2,'C':.3,'G':.3,'T':.2}
    kappa = 2
    rate_matrix_object = RateMatrix.get_unscaled_hky85_rate_matrix(distribution, kappa)
    rate_matrix_object.normalize()
    rate_matrix = rate_matrix_object.get_dictionary_rate_matrix()
    path_length = 2
    initial_state = 'A'
    terminal_state = 'C'
    states = 'ACGT'
    uniformization_events = get_uniformization_sample(initial_state, terminal_state, states, path_length, rate_matrix)
    print uniformization_events

def main():
    #demo_rejection_sampling()
    demo_discrete_path_sampling()
    #demo_uniformization()
    

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    #parser.add_option('-v', '--verbose', action='store_true', dest='verbose', default=False)
    #parser.add_option('-o', '--output', dest='output_filename', metavar='FILE', help='output file')
    parser.add_option('--test', action='store_true', dest='test', default=False)
    options, args = parser.parse_args()
    if options.test:
        suite = unittest.TestLoader().loadTestsFromTestCase(TestPathSampler)
        unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        main()


