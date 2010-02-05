"""
Convert a rate matrix to a transition matrix.
Use matrix exponentiation to get a transition matrix
from a rate matrix and a time or branch length.
"""

# standard modules
from StringIO import StringIO
import unittest
import math

# extension modules
import numpy as np
import scipy

# nonstandard modules
import CodonFrequency
import Codon
import MatrixUtil
import Util


class RateMatrixError(Exception):
    pass


def almost_equal(a, b, eps=0.000000001):
    return abs(a - b) < eps

def get_sample_codon_rate_matrix():
    """
    Get the Goldman-Yang 1994 codon rate matrix used by Eric Stone and Asger Hobolth for path simulation.
    """
    cf = CodonFrequency.codon_frequency_b
    distribution = dict((codon, cf.codon_to_non_stop_proportion(codon)) for codon in Codon.g_non_stop_codons)
    kappa = 2
    omega = .01
    codon_rate_matrix = get_gy94_rate_matrix(distribution, kappa, omega)
    return codon_rate_matrix

def get_stationary_distribution(row_major_rate_matrix):
    """
    @return: a list of stationary frequencies
    """
    # convert the transition matrix to numpy array form
    numpy_array = np.array(row_major_rate_matrix)
    # get the eigensystem
    left = True
    right = False
    w, vl = scipy.linalg.eig(numpy_array, None, left, right)
    # find the eigenvalue with the smallest absolute value
    best_eigenvalue_size, best_index = min((abs(value), i) for i, value in enumerate(w))
    # get the stationary distribution
    scaled_stationary_distribution = list(vl.T[best_index])
    scaled_total = sum(scaled_stationary_distribution)
    stationary_distribution = [value / scaled_total for value in scaled_stationary_distribution]
    # assert the validity of the stationary distribution
    if min(stationary_distribution) < 0:
        raise RateMatrixError('found a negative element of the stationary eigenvector: %s' % min(stationary_distribution))
    return stationary_distribution

def get_eigendecomposition(row_major_matrix):
    """
    @return: a numpy matrix of right eigenvectors with correct handling of repeated eigenvalues
    """
    p = get_stationary_distribution(row_major_matrix)
    prematrix = scipy.diag([math.sqrt(x) for x in p])
    postmatrix = scipy.diag([1/math.sqrt(x) for x in p])
    M = np.array(row_major_matrix)
    # define a symmetric matrix
    S = np.dot(np.dot(prematrix, M), postmatrix)
    # get the spectral decomposition of the symmetric matrix
    w, vr = scipy.linalg.eigh(S)
    #print 'this should be the identity matrix:'
    #print np.dot(vr, vr.T)
    # reconstruct the eigenvector matrix of the original matrix
    U = np.dot(postmatrix, vr)
    #print 'this should also be the identity matrix:'
    #U_inv = np.dot(vr.T, prematrix)
    #print np.dot(U, U_inv)
    return w, U

def get_likelihood(tree, rate_matrix_object):
    """
    Get the likelihood given an aligned column and a tree.
    @param tree: a tree with branch lengths and leaf states
    @param rate_matrix_object: a rate matrix object
    @return: a likelihood
    """
    rate_matrix_object.validate_states(set(node.state for node in tree.gen_tips()))
    add_probabilities(tree, rate_matrix_object)
    likelihood = 0
    for state, p in zip(rate_matrix_object.states, rate_matrix_object.get_stationary_distribution()):
        likelihood += p * tree.root.state_to_subtree_prob[state]
    return likelihood

def add_probabilities(tree, rate_matrix_object):
    """
    Augment each node in the tree with conditional probabilities.
    This is used for calculating likelihood and for simulating ancestral states.
    @param tree: a tree with branch lengths and leaf states
    @param rate_matrix_object: a rate matrix object with states and dictionary transition matrices
    """
    rate_matrix_object.validate_states(set(node.state for node in tree.gen_tips()))
    for node in tree.postorder():
        node.state_to_subtree_prob = {}
        for state in rate_matrix_object.states:
            if not node.children:
                if state == node.state:
                    node.state_to_subtree_prob[state] = 1
                else:
                    node.state_to_subtree_prob[state] = 0
            else:
                probability = 1
                for child in node.children:
                    subtree_probability = 0
                    transition_matrix = rate_matrix_object.get_dictionary_transition_matrix(child.get_branch_length())
                    if child.children:
                        for child_state in rate_matrix_object.states:
                            conditional_transition_probability = transition_matrix[(state, child_state)]
                            conditional_subtree_probability = child.state_to_subtree_prob[child_state]
                            subtree_probability += conditional_transition_probability * conditional_subtree_probability
                    else:
                        subtree_probability = transition_matrix[(state, child.state)]
                    probability *= subtree_probability
                node.state_to_subtree_prob[state] = probability


class FastRateMatrix:
    """
    This is an experimental variant of RateMatrix with a reduced interface.
    """

    def __init__(self, row_major_matrix, ordered_states):
        """
        @param row_major_matrix: the rate matrix in row major format
        @param ordered_states the names of the the states
        """
        # store the state names
        self.states = ordered_states
        # get the stationary distribution
        self.stationary_distribution = get_stationary_distribution(row_major_matrix)
        p = self.stationary_distribution
        # get the expected rate of the rate matrix
        self.rate = -sum(x*row_major_matrix[i][i] for i, x in enumerate(self.stationary_distribution))
        # get the numpy matrix
        self.matrix = np.array(row_major_matrix)
        M = self.matrix
        # get the diagonal matrix representing the square root of the stationary frequencies
        D = scipy.diag([math.sqrt(x) for x in p])
        # get the inverse of the diagonal matrix
        D_inv = scipy.diag([1/math.sqrt(x) for x in p])
        # define a symmetric matrix
        S = np.dot(np.dot(D, M), D_inv)
        # get the spectral decomposition of the symmetrice matrix
        w, vr = scipy.linalg.eigh(S)
        # reconstruct the eigenvector matrix of the original matrix
        self.diagonal = [x/self.rate for x in w]
        self.pre_matrix = np.dot(D_inv, vr)
        self.post_matrix = np.dot(vr.T, D)
        # define the cache
        self.distance_to_transition_matrix = {}

    def get_expected_times(self, initial_state, final_state, t):
        """
        @param initial_state: the state at one end of the evolutionary branch
        @param final_state: the state at the other end of the evolutionary branch
        @param t: the amount of time between the initial and final states
        @return: a list that defines the expected time spent in each state
        """
        # define the eigenvector matrices
        #eigenvalues, U = get_eigendecomposition(self.get_row_major_rate_matrix())
        #U_inv = scipy.linalg.inv(U)
        """
        # validate the eigendecomposition
        D = scipy.diag(eigenvalues)
        lhs = np.dot(self.matrix, U)
        rhs = np.dot(U, D)
        assert np.all(lhs == rhs), (lhs, rhs)
        """
        # define a diagonal matrix using exponentiated eigenvalues
        eigenvalues = [x * self.rate for x in self.diagonal]
        exponentiated_eigenvalues = [math.exp(t*x) for x in eigenvalues]
        exponentiated_diagonal = scipy.diag(exponentiated_eigenvalues)
        # define the transition matrix for checking the calculation
        transition_matrix = np.dot(np.dot(self.pre_matrix, exponentiated_diagonal), self.post_matrix)
        U = self.pre_matrix
        U_inv = self.post_matrix
        """
        # validate the transition matrix
        lhs = transition_matrix
        rhs = scipy.linalg.expm(self.matrix*t)
        assert np.all(lhs == rhs), (lhs, rhs)
        """
        # get the unscaled wait time for each state
        matrix_size = len(exponentiated_eigenvalues)
        unscaled_wait_times = []
        for wait_state in range(matrix_size):
            # calculate the integral
            total = 0
            for i in range(matrix_size):
                partial_value = U[initial_state][i] * U_inv[i][wait_state]
                subtotal = 0
                for j in range(matrix_size):
                    value = partial_value * U[wait_state][j] * U_inv[j][final_state]
                    if almost_equal(eigenvalues[i], eigenvalues[j]):
                        value *= t * exponentiated_eigenvalues[i]
                    else:
                        dy = exponentiated_eigenvalues[i] - exponentiated_eigenvalues[j]
                        dx = eigenvalues[i] - eigenvalues[j]
                        value *= (dy / dx)
                    subtotal += value
                total += subtotal
            unscaled_wait_times.append(float(total))
        # scale the wait times
        scaling_factor = 1 / transition_matrix[initial_state][final_state]
        scaled_wait_times = [x * scaling_factor for x in unscaled_wait_times]
        # we should be able to check the calculation using the transition matrix
        total_wait_time = sum(scaled_wait_times)
        assert almost_equal(t, total_wait_time), (t, total_wait_time)
        return scaled_wait_times

    def get_row_major_rate_matrix(self):
        diagonal_matrix = scipy.diag(self.diagonal)
        numpy_matrix = np.dot(np.dot(self.pre_matrix, diagonal_matrix), self.post_matrix)
        row_major_rate_matrix = []
        m = len(self.states)
        for i in range(m):
            row = [self.matrix[i][j] for j in range(m)]
            row_major_rate_matrix.append(row)
        return row_major_rate_matrix

    def create_dictionary_transition_matrix(self, t):
        """
        Create a transition matrix given a time between initial and final states.
        @param t: the amount of time between the initial and final states
        @return: a dictionary transition matrix
        """
        distance = self.rate * t
        diagonal_matrix = scipy.diag([math.exp(distance*x) for x in self.diagonal])
        numpy_transition_matrix = np.dot(np.dot(self.pre_matrix, diagonal_matrix), self.post_matrix)
        d = {}
        for i, a in enumerate(self.states):
            for j, b in enumerate(self.states):
                d[(a,b)] = numpy_transition_matrix[i][j]
        return d

    def get_dictionary_transition_matrix(self, t):
        """
        This is a memoized version of create_dictionary_transition_matrix().
        @param t: the amount of time between the initial and final states
        @return: a dictionary transition matrix
        """
        distance = self.rate * t
        if distance not in self.distance_to_transition_matrix:
            self.distance_to_transition_matrix[distance] = self.create_dictionary_transition_matrix(t)
        return self.distance_to_transition_matrix[distance]

    def validate_states(self, states):
        """
        Verify that each of the input states is a valid state.
        @param states: states that should be valid
        @raise RateMatrixError: a state is invalid
        """
        bad_states = set(states) - set(self.states)
        if bad_states:
            example_bad_state = bad_states.pop()
            raise RateMatrixError('invalid state: %s' % example_bad_state)

    def normalize(self):
        self.rate = 1.0

    def set_rate(self, rate):
        self.rate = rate

    def get_stationary_distribution(self):
        return self.stationary_distribution




class RateMatrix:
    """
    Cache transition matrices for various branch lengths.
    """

    def __init__(self, row_major_matrix, ordered_states):
        """
        @param row_major_matrix: the rate matrix in row major form
        @param ordered_states: the state labels in an order that corresponds to the row_major_matrix order
        """
        # do some sanity checks
        assert row_major_matrix
        nrows = len(row_major_matrix)
        ncols = len(row_major_matrix[0])
        assert nrows == ncols
        assert len(ordered_states) == nrows
        assert len(ordered_states) == ncols
        # create the dictionary rate matrix
        self.dictionary_rate_matrix = MatrixUtil.row_major_to_dict(row_major_matrix, ordered_states, ordered_states)
        # create the numpy matrix
        self.numpy_rate_matrix = np.array(row_major_matrix)
        # cache some stuff
        self.row_major_rate_matrix = row_major_matrix
        self.states = ordered_states
        self.state_to_index = dict((state, i) for i, state in enumerate(self.states))
        # create some more or less uninitialized variables
        self.branch_length_to_numpy_transition_matrix = {}
        self.branch_length_to_dictionary_transition_matrix = {}
        self.stationary_distribution = None

    def get_row_major_rate_matrix(self):
        return self.row_major_rate_matrix

    def get_dictionary_rate_matrix(self):
        return self.dictionary_rate_matrix

    def rescale(self, scaling_factor):
        """
        Rescale the rate matrix.
        Some cached values are invalidated.
        The stationary distribution is not affected.
        @param scaling_factor: each element of the rate matrix is multiplied by this factor
        """
        # TODO the object should keep track of the current scaling,
        # and this function should just change that number;
        # the cached transition matrices could be preserved.
        #
        # handle a degenerate case
        if scaling_factor == 1:
            return
        # invalidate all cached transition matrices
        self.branch_length_to_numpy_transition_matrix = {}
        self.branch_length_to_dictionary_transition_matrix = {}
        # modify the row major rate matrix
        self.row_major_rate_matrix = [[x * scaling_factor for x in row] for row in self.row_major_rate_matrix]
        # regenerate the rate matrices using different formats
        self.dictionary_rate_matrix = MatrixUtil.row_major_to_dict(self.row_major_rate_matrix, self.states, self.states)
        self.numpy_rate_matrix = np.array(self.row_major_rate_matrix)

    def normalize(self):
        """
        Rescale the rate matrix to have an expected rate of one.
        """
        scaling_factor = 1.0 / self.get_expected_rate()
        self.rescale(scaling_factor)
        eps = .000001
        if abs(self.get_expected_rate() - 1.0) > eps:
            raise RateMatrixError('the expected rate should be 1.0 after normalization')

    def get_stationary_distribution(self):
        if self.stationary_distribution is None:
            self.stationary_distribution = get_stationary_distribution(self.row_major_rate_matrix)
        return self.stationary_distribution

    def get_numpy_transition_matrix(self, t):
        """
        @param t: the time or distance over which the transition occurs
        @return: a numpy transition matrix
        """
        transition_matrix = self.branch_length_to_numpy_transition_matrix.get(t, None)
        if transition_matrix is not None:
            return transition_matrix
        transition_matrix = scipy.linalg.expm(self.numpy_rate_matrix * t)
        self.branch_length_to_numpy_transition_matrix[t] = transition_matrix
        return transition_matrix

    def _create_dictionary_transition_matrix(self, t):
        """
        @param t: the time or distance over which the transition occurs
        @return: a dictionary transition matrix
        """
        # create the dictionary from the numpy transition matrix
        numpy_transition_matrix = self.get_numpy_transition_matrix(t)
        rows, cols = numpy_transition_matrix.shape
        assert rows == cols
        assert rows == len(self.states)
        dictionary_transition_matrix = {}
        for row_index in range(rows):
            a = self.states[row_index]
            for col_index in range(cols):
                b = self.states[col_index]
                dictionary_transition_matrix[(a, b)] = numpy_transition_matrix[row_index, col_index]
        # add the dictionary to the cache and return the dictionary
        self.branch_length_to_dictionary_transition_matrix[t] = dictionary_transition_matrix
        return dictionary_transition_matrix

    def get_dictionary_transition_matrix(self, t):
        """
        @param t: the time or distance over which the transition occurs
        @return: a dictionary transition matrix
        """
        # see if we have the dictionary cached already
        dictionary_transition_matrix = self.branch_length_to_dictionary_transition_matrix.get(t, None)
        if dictionary_transition_matrix is None:
            dictionary_transition_matrix = self._create_dictionary_transition_matrix(t)
            self.branch_length_to_dictionary_transition_matrix[t] = dictionary_transition_matrix
        return dictionary_transition_matrix

    def validate_states(self, states):
        """
        Verify that each of the input states is a valid state.
        @param states: states that should be valid
        @raise RateMatrixError: a state is invalid
        """
        bad_states = set(states) - set(self.states)
        if bad_states:
            example_bad_state = bad_states.pop()
            raise RateMatrixError('invalid state: %s' % example_bad_state)

    def add_probabilities(self, tree):
        """
        Augment each node in the tree with conditional probabilities.
        This is used for calculating likelihood and for simulating ancestral states.
        @param tree: a tree with branch lengths and leaf states
        """
        #FIXME obsolete in favor of a generic function
        self.validate_states(set(node.state for node in tree.gen_tips()))
        for node in tree.postorder():
            node.state_to_subtree_prob = {}
            for state in self.states:
                if not node.has_children():
                    if state == node.state:
                        node.state_to_subtree_prob[state] = 1
                    else:
                        node.state_to_subtree_prob[state] = 0
                else:
                    probability = 1
                    for child in node.gen_children():
                        subtree_probability = 0
                        transition_matrix = self.get_dictionary_transition_matrix(child.get_branch_length())
                        if child.has_children():
                            for child_state in self.states:
                                conditional_transition_probability = transition_matrix[(state, child_state)]
                                conditional_subtree_probability = child.state_to_subtree_prob[child_state]
                                subtree_probability += conditional_transition_probability * conditional_subtree_probability
                        else:
                            subtree_probability = transition_matrix[(state, child.state)]
                        probability *= subtree_probability
                    node.state_to_subtree_prob[state] = probability

    def simulate_states(self, tree):
        """
        Augment each node in the tree with a simulated state.
        @param tree: a tree with branch lengths
        """
        for node in tree.preorder():
            if node is tree.root:
                weights = self.get_stationary_distribution()
            else:
                d = self.get_dictionary_transition_matrix(node.get_branch_length())
                parent_state = node.get_parent().state
                weights = [d[(parent_state, state)] for state in self.states]
            node.state = Util.weighted_choice(zip(weights, self.states))

    def simulate_ancestral_states(self, tree):
        """
        Augment each internal node in the tree with a simulated state.
        This simulation is conditional on the states at the leaves.
        @param tree: a tree with branch lengths and leaf states
        """
        self.validate_states(set(node.state for node in tree.gen_tips()))
        self.add_probabilities(tree)
        for node in tree.gen_internal_nodes_preorder():
            weights = []
            for state, p in zip(self.states, self.get_stationary_distribution()):
                conditional_subtree_probability = node.state_to_subtree_prob[state]
                if node is tree.root:
                    conditional_transition_probability = p
                else:
                    parent_state = node.parent.state
                    d = self.get_dictionary_transition_matrix(node.get_branch_length())
                    conditional_transition_probability = d[(parent_state, state)]
                weights.append(conditional_subtree_probability*conditional_transition_probability)
            node.state = Util.weighted_choice(zip(weights, self.states))

    def get_likelihood(self, tree):
        """
        Get the likelihood given an aligned column and a tree.
        @param tree: a tree with branch lengths and leaf states
        @return: a likelihood
        """
        #FIXME obsolete in favor of a generic function
        self.validate_states(set(node.state for node in tree.gen_tips()))
        self.add_probabilities(tree)
        likelihood = 0
        for state, p in zip(self.states, self.get_stationary_distribution()):
            likelihood += p * tree.root.state_to_subtree_prob[state]
        return likelihood

    def get_expected_rate(self):
        """
        Get the expected substitution rate.
        @return: the expected substitution rate
        """
        expected_rate = 0
        for state, p in zip(self.states, self.get_stationary_distribution()):
            rate = -self.dictionary_rate_matrix[(state, state)]
            expected_rate += p * rate
        return expected_rate


def get_unscaled_hky85_rate_matrix(distribution, kappa):
    """
    @param distribution: a dictionary mapping a nucleotide to its stationary frequency
    @param kappa: the transition / transversion substitution rate ratio
    @return: a nucleotide rate matrix object
    """
    assert len(distribution) == 4
    assert set(distribution) == set('ACGT')
    assert abs(sum(distribution.values()) - 1.0) < .0000001
    # Create the off-diagonal elements of the unscaled rate matrix.
    rate_matrix = {}
    for na in distribution:
        for nb, probability in distribution.items():
            if na != nb:
                rate = probability
                if na+nb in ('AG', 'GA', 'CT', 'TC'):
                    rate *= kappa
                rate_matrix[(na, nb)] = rate
    # Create the diagonal elements such that each row in the rate matrix sums to zero.
    for na in distribution:
        rate = sum(rate_matrix[(na, nb)] for nb in distribution if nb != na)
        rate_matrix[(na, na)] = -rate
    # Convert the dictionary rate matrix to a row major rate matrix
    ordered_states = list('ACGT')
    row_major_rate_matrix = MatrixUtil.dict_to_row_major(rate_matrix, ordered_states, ordered_states)
    rate_matrix_object = RateMatrix(row_major_rate_matrix, ordered_states)
    return rate_matrix_object

def get_gy94_rate_matrix(distribution, kappa, omega):
    """
    This codon rate matrix was described by Goldman and Yang in 1994.
    @param distribution: a dictionary mapping a codon to its stationary frequency
    @param kappa: the transition / transversion substitution rate ratio
    @param omega: the non-synonymous / synonymous substitution rate ratio
    @return: a codon rate matrix in convenient dictionary form
    """
    # TODO make this return a rate matrix object instead of a dictionary
    # There should be 61 codons because stop codons are not included.
    assert len(distribution) == 61
    # The sum of the stationary probabilities should be one.
    assert abs(sum(distribution.values()) - 1.0) < .0000001
    # The codons should be valid.
    assert set(distribution) == set(Codon.g_non_stop_codons)
    # Create the off-diagonal elements of the unscaled rate matrix.
    rate_matrix = {}
    for ca in distribution:
        for cb, probability in distribution.items():
            if ca != cb:
                transition_count = 0
                transversion_count = 0
                for a, b in zip(ca, cb):
                    if a != b:
                        if a+b in ('AG', 'GA', 'CT', 'TC'):
                            transition_count += 1
                        else:
                            transversion_count += 1
                rate = 0
                if transition_count + transversion_count == 1:
                    rate = probability
                    if transition_count:
                        rate *= kappa
                    if Codon.g_codon_to_aa_letter[ca] != Codon.g_codon_to_aa_letter[cb]:
                        rate *= omega
                rate_matrix[(ca, cb)] = rate
    # Create the diagonal elements such that each row in the rate matrix sums to zero.
    for ca in distribution:
        rate = sum(rate_matrix[(ca, cb)] for cb in distribution if cb != ca)
        rate_matrix[(ca, ca)] = -rate
    # Scale the rate matrix so that t substitutions are expected in t time units.
    mu = -sum(rate_matrix[(codon, codon)] * probability for codon, probability in distribution.items())
    for ca in distribution:
        for cb in distribution:
            rate_matrix[(ca, cb)] /= mu
    return rate_matrix


def get_jukes_cantor_rate_matrix():
    states = 'ACGT'
    rate_matrix = {}
    for src in states:
        for dst in states:
            if src != dst:
                rate_matrix[(src, dst)] = 1.0 / 3.0
    for src in states:
        rate_matrix[(src, src)] =  -sum(rate_matrix[(src, dst)] for dst in states if src != dst)
    return rate_matrix

def get_jukes_cantor_transition_matrix(t):
    states = 'ACGT'
    transition_matrix = {}
    for src in states:
        for dst in states:
            if src == dst:
                p = (1.0/4.0) + (3.0/4.0)*math.exp(-(4.0/3.0)*t)
            else:
                p = (1.0/4.0) - (1.0/4.0)*math.exp(-(4.0/3.0)*t)
            transition_matrix[(src, dst)] = p
    return transition_matrix


class TestRateMatrix(unittest.TestCase):

    def test_rate_matrix(self):
        """
        Compare the transition matrix found using a RateMatrix object to the transition matrix found analytically.
        """
        # get the row major Jukes-Cantor rate matrix
        dictionary_rate_matrix = get_jukes_cantor_rate_matrix()
        ordered_states = list('ACGT')
        row_major_rate_matrix = MatrixUtil.dict_to_row_major(dictionary_rate_matrix, ordered_states, ordered_states)
        # get the rate matrix object
        rate_matrix_object = RateMatrix(row_major_rate_matrix, ordered_states)
        # get the transition matrix
        t = 2.0
        calculated_transition_matrix = rate_matrix_object.get_dictionary_transition_matrix(t)
        analytical_transition_matrix = get_jukes_cantor_transition_matrix(t)
        self.assertEqual(len(calculated_transition_matrix), 16)
        self.assertEqual(len(analytical_transition_matrix), 16)
        for key, value in calculated_transition_matrix.items():
            self.assertAlmostEqual(analytical_transition_matrix[key], value)

    def test_arbitrary_reversible_asymmetric_rate_matrices(self):
        """
        Compare the transition matrices found using different numerical methods.
        """
        # create an arbitrary but asymmetric and reversible rate matrix
        row_major_rate_matrix = [
                [-2, 1, 1],
                [2, -3, 1],
                [2, 1, -3]]
        states = ['a', 'b', 'c']
        # define the amount of time over which the transition happens
        t = 2.0
        # get the dictionary transition matrix from the RateMatrix object
        rate_matrix_object = RateMatrix(row_major_rate_matrix, states)
        d1 = rate_matrix_object.get_dictionary_transition_matrix(t)
        # get the dictionary transition matrix from the FastRateMatrix object
        fast_rate_matrix_object = FastRateMatrix(row_major_rate_matrix, states)
        d2 = fast_rate_matrix_object.get_dictionary_transition_matrix(t)
        # compare the elements of the dictionaries
        transitions = [(a, b) for a in states for b in states]
        for transition in transitions:
            self.assertAlmostEqual(d1[transition], d2[transition])

    def test_stationary_distribution_a(self):
        """
        Test the stationary distribution found by a RateMatrix object.
        """
        # create an arbitrary but asymmetric and reversible rate matrix
        row_major_rate_matrix = [
                [-2, 1, 1],
                [2, -3, 1],
                [2, 1, -3]]
        # define the actual stationary distribution
        actual_stationary_distribution = [.5, .25, .25]
        # get the calculated stationary distribution
        observed_stationary_distribution = get_stationary_distribution(row_major_rate_matrix)
        # corresponding elements should be almost equal
        for a, b in zip(actual_stationary_distribution, observed_stationary_distribution):
            self.assertAlmostEqual(a, b)


def codon_rate_matrix_to_html_string(rate_matrix):
    codons = list(sorted(Codon.g_non_stop_codons))
    arr = []
    first_row_arr = ['<tr><td></td>'] + ['<td>%s</td>' % codon for codon in codons] + ['</tr>']
    arr.append(''.join(first_row_arr))
    for ca in Codon.g_non_stop_codons:
        row_arr = ['<tr><td>%s</td>' % ca] + ['<td>%.3f</td>' % rate_matrix[(ca, cb)] for cb in codons] + ['</tr>']
        arr.append(''.join(row_arr))
    return '\n'.join(arr)

def demo_codon_table():
    import CodonFrequency
    print '<html>'
    print '<head><style type="text/css">td{font-size:x-small;}</style></head>'
    print '<body>'
    for cf in (CodonFrequency.codon_frequency_a, CodonFrequency.codon_frequency_b):
        distribution = dict((codon, cf.codon_to_non_stop_proportion(codon)) for codon in Codon.g_non_stop_codons)
        rate_matrix = get_gy94_rate_matrix(distribution, 2, .01)
        table_guts = codon_rate_matrix_to_html_string(rate_matrix)
        print '<table>' + table_guts + '</table><br/><br/>'
    print '</body>'
    print '</html>'

def demo_rate_matrix_eigendecomposition():
    # get a sample codon rate matrix that is in the form of a dictionary
    rate_matrix = get_sample_codon_rate_matrix()
    states = list(sorted(Codon.g_non_stop_codons))
    # convert the rate matrix to row major form
    row_major_matrix = MatrixUtil.dict_to_row_major(rate_matrix, states, states)
    # show the eigensystem of the row major matrix
    print 'row sums:'
    for row in row_major_matrix:
        print sum(row)
    # convert the transition matrix to numpy array form
    numpy_array = np.array(row_major_matrix)
    # get the eigensystem
    left = True
    right = True
    w, vl, vr = scipy.linalg.eig(numpy_array, None, left, right)
    # show the eigensystem
    print 'eigenvalues:'
    print w
    print 'left eigenvectors:'
    print vl.T
    print 'right eigenvectors:'
    print vr.T


def demo_rate_matrix_stationary_distribution():
    # get a sample codon rate matrix that is in the form of a dictionary
    rate_matrix = get_sample_codon_rate_matrix()
    states = list(sorted(Codon.g_non_stop_codons))
    # convert the rate matrix to row major form
    row_major_matrix = MatrixUtil.dict_to_row_major(rate_matrix, states, states)
    # show the eigensystem of the row major matrix
    print 'row sums:'
    for row in row_major_matrix:
        print sum(row)
    # convert the transition matrix to numpy array form
    numpy_array = np.array(row_major_matrix)
    # get the eigensystem
    left = True
    right = False
    w, vl = scipy.linalg.eig(numpy_array, None, left, right)
    # find the eigenvalue with the smallest absolute value
    best_eigenvalue_size, best_index = min((abs(value), i) for i, value in enumerate(w))
    print 'smallest absolute value of an eigenvalue:', best_eigenvalue_size
    # get the stationary distribution
    scaled_stationary_distribution = list(vl.T[best_index])
    scaled_total = sum(scaled_stationary_distribution)
    stationary_distribution = [value / scaled_total for value in scaled_stationary_distribution]
    # assert the validity of the stationary distribution
    if min(stationary_distribution) < 0:
        raise RateMatrixError('found a negative element of the stationary eigenvector: %s' % min(stationary_distribution))
    print 'stationary distribution:'
    for value in stationary_distribution:
        print value

def demo_eigenvector_orthonormality():
    """
    Demo code to find the expected state distribution on a path.
    """
    """
    row_major = [
            [-3, 1, 1, 1],
            [1, -3, 1, 1],
            [1, 1, -3, 1],
            [1, 1, 1, -3]]
    """
    row_major = [
            [-.3, .1, .1, .1],
            [.2, -.4, .1, .1],
            [.2, .1, -.4, .1],
            [.2, .1, .1, -.4]]
    M = np.array(row_major)
    print 'rate matrix (M):'
    print M
    matrix_size = len(row_major)
    print 'M is', matrix_size, 'x', matrix_size
    eigenvalues, U = get_eigendecomposition(row_major)
    D = scipy.diag(eigenvalues)
    print len(set(eigenvalues)), 'distinct eigenvalues'
    print 'diagonal eigenvalue matrix (D):'
    print D
    print 'right eigenvector matrix (U):'
    print U
    print 'transpose of U:'
    print U.T
    U_inv = scipy.linalg.inv(U)
    print 'inverse of U:'
    print U_inv
    print 'M x U:'
    print np.dot(M, U)
    print 'U x D:'
    print np.dot(U, D)
    print 'U x U^T:'
    print np.dot(U, U.T)
    exp_eigenvalues = [math.exp(x) for x in eigenvalues]
    exp_D = scipy.diag(exp_eigenvalues)
    transition_matrix = np.dot(np.dot(U, exp_D), U_inv)
    print 'transition matrix:'
    print transition_matrix
    print 'reconstruction from the eigendecomposition:'
    print np.dot(np.dot(U, D), U_inv)
    initial_state = 0
    final_state = 1
    unscaled_wait_times = []
    for wait_state in range(matrix_size):
        # calculate the integral
        total = 0
        for i in range(matrix_size):
            for j in range(matrix_size):
                value = U[initial_state][i] * U_inv[i][wait_state] * U[wait_state][j] * U_inv[j][final_state]
                if almost_equal(eigenvalues[i], eigenvalues[j]):
                    value *= exp_eigenvalues[i]
                else:
                    dy = exp_eigenvalues[i] - exp_eigenvalues[j]
                    dx = eigenvalues[i] - eigenvalues[j]
                    value *= (dy / dx)
                total += value
        unscaled_wait_times.append(float(total))
    # scale the wait times
    scaling_factor = 1 / transition_matrix[initial_state][final_state]
    scaled_wait_times = [x * scaling_factor for x in unscaled_wait_times]
    # we should be able to check the calculation using the transition matrix
    total_wait_time = sum(scaled_wait_times)
    print 'total wait time:', total_wait_time
    print 'expected total wait time:', 1
    print 'scaled wait times:'
    print scaled_wait_times

def main():
    #demo_rate_matrix_stationary_distribution()
    #demo_rate_matrix_eigendecomposition()
    demo_eigenvector_orthonormality()

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    #parser.add_option('-v', '--verbose', action='store_true', dest='verbose', default=False)
    #parser.add_option('-o', '--output', dest='output_filename', metavar='FILE', help='output file')
    parser.add_option('--test', action='store_true', dest='test', default=False)
    options, args = parser.parse_args()
    if options.test:
        suite = unittest.TestLoader().loadTestsFromTestCase(TestRateMatrix)
        unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        main()


