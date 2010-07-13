"""
Do things with transition matrices.
"""

import unittest

import numpy as np

import Util
import DiscreteEndpoint


class TransitionObject:
    def get_transition_probability(self, source, sink, distance):
        """
        @param source: the source state index
        @param sink: the sink state index
        @param distance: the number of steps
        @return: transition probability
        """
        raise NotImplementedError()
    def get_stationary_probability(self, state):
        raise NotImplementedError()
    def get_stationary_distribution(self):
        raise NotImplementedError()
    def get_ntransitions_expected(self, source, sink, distance):
        raise NotImplementedError()
    def get_nstates(self):
        raise NotImplementedError()


class MatrixTransitionObject:
    """
    This is like a transition matrix.
    Matrix powers could use some caching eventually.
    """

    def __init__(self, T):
        """
        @param T: a right stochastic matrix as a numpy array
        """
        self.T = T
        self.stationary_distribution = get_stationary_distribution(self.T)

    def get_nstates(self):
        return len(self.stationary_distribution)

    def get_transition_probability(self, source, sink, distance=1):
        """
        @param source: the source state index
        @param sink: the sink state index
        @param distance: the number of steps
        @return: transition probability
        """
        if distance < 1:
            raise ValueError('expected a positive integer')
        return np.linalg.matrix_power(self.T, distance)[source, sink]

    def get_stationary_probability(self, state):
        return self.stationary_distribution[state]

    def get_stationary_distribution(self):
        """
        @return: a stochastic vector as a list
        """
        return self.stationary_distribution

    def get_ntransitions_expected(self, source, sink, distance):
        raise NotImplementedError()


class UniformTransitionObject:
    """
    This is like a transition matrix.
    Each state is assumed to have the same stationary probability.
    Randomization occurs between discrete steps according to a poisson process.
    """

    def __init__(self, prandom, nstates, cache_size=0):
        """
        @param prandom: the probability of randomization per step
        @param nstates: the number of possible states
        @param cache_size: save this many observations
        """
        if not (0 <= prandom <= 1):
            raise ValueError('expected a probability')
        self.prandom = prandom
        self.nstates = nstates
        self._get_ntrans = Util.Cache(self._get_uncached_ntrans, cache_size)

    def get_nstates(self):
        return self.nstates

    def get_transition_probability(self, source, sink, distance=1):
        """
        @param source: the source state index
        @param sink: the sink state index
        @param distance: the number of steps
        @return: transition probability
        """
        prandom_total = 1 - (1 - self.prandom)**distance
        if source == sink:
            return prandom_total / self.nstates + (1 - prandom_total)
        else:
            return prandom_total / self.nstates

    def get_stationary_distribution(self):
        """
        @return: a stochastic vector as a list
        """
        return [1.0 / self.nstates] * self.nstates

    def get_stationary_probability(self, state):
        return 1.0 / self.nstates

    def get_ntransitions_expected(self, source, sink, distance):
        e_same, e_different = self._get_ntrans(distance)
        if source == sink:
            return e_same
        else:
            return e_different

    def _get_uncached_ntrans(self, distance):
        return DiscreteEndpoint.get_expected_transitions_binomial(self.prandom, self.nstates, distance)


def get_stationary_distribution(transition_matrix):
    """
    @param transition_matrix: a right stochastic matrix like a numpy array
    @return: a stochastic vector as a list
    """
    T = np.array(transition_matrix)
    # Do validation
    nrows, ncols = T.shape
    if nrows != ncols:
        raise ValueError('expected a square transition matrix')
    if not np.allclose(np.sum(T, 1), np.ones(ncols)):
        raise ValueError('expected a right stochastic transition matrix')
    # We want a left eigenvector of T.
    # Numpy's eig gives only the right eigenvectors,
    # so use the transpose of T.
    w, VT = np.linalg.eig(T.T)
    # We want the eigenvector that corresponds to the eigenvalue of 1.
    # No eigenvalue should be greater than 1.
    best_w, best_v = max((x, v) for x, v in zip(w, VT.T))
    # The eigenvector might have tiny imaginary parts, so remove them.
    best_v = np.array([abs(x) for x in best_v])
    # Force the elements of the dominant eigenvector to sum to one.
    best_v /= np.sum(best_v)
    return best_v.tolist()

def get_uniform_transition_matrix(prandom, nstates, ntransitions=1):
    """
    @param prandom: probability of randomization per transition
    @param nstates: the number of states
    @param ntransitions: the transition matrix is over this many transitions
    @return: a numpy array representing the transition matrix
    """
    if not (0 <= prandom <= 1):
        raise ValueError('expected a probability')
    if ntransitions < 1:
        raise ValueError('expected a positive integer')
    prandom_total = 1 - (1 - prandom)**ntransitions
    T = np.ones((nstates, nstates)) * (prandom_total/nstates)
    T += np.diag([1-prandom_total]*nstates)
    return T


class TestTransitionMatrix(unittest.TestCase):

    def test_stationary_distribution_a(self):
        """
        Use the Wikipedia example from el gato, el reloj y el raton.
        """
        T = np.array([
                [0.00, 0.00, 0.50, 0.00, 0.50],
                [0.00, 0.00, 1.00, 0.00, 0.00],
                [0.25, 0.25, 0.00, 0.25, 0.25],
                [0.00, 0.00, 0.50, 0.00, 0.50],
                [0.00, 0.00, 0.00, 0.00, 1.00]])
        observed = get_stationary_distribution(T)
        expected = [0, 0, 0, 0, 1]
        self.assertTrue(np.allclose(observed, expected))

    def test_stationary_distribution_b(self):
        """
        Use an example from Applied Linear Algebra by Noble and Daniel.
        """
        T = np.array([
            [0.8, 0.1, 0.1],
            [0.2, 0.7, 0.1],
            [0.1, 0.3, 0.6]])
        observed = get_stationary_distribution(T)
        expected = [0.45, 0.35, 0.20]
        self.assertTrue(np.allclose(observed, expected), str(observed))

    def test_stationary_distribution_c(self):
        """
        Use an irreversible rock-scissors-paper-like matrix.
        This matrix has imaginary eigenvalues.
        """
        T = np.array([
            [0.1, 0.8, 0.1],
            [0.1, 0.1, 0.8],
            [0.8, 0.1, 0.1]])
        observed = get_stationary_distribution(T)
        expected = np.ones(3)/3.0
        self.assertTrue(np.allclose(observed, expected), str(observed))

    def test_stationary_distribution_d(self):
        """
        Compare stationary distributions provided by the objects.
        """
        prandom = .1
        nstates = 4
        M = np.array([
            [.925, .025, .025, .025],
            [.025, .925, .025, .025],
            [.025, .025, .925, .025],
            [.025, .025, .025, .925]])
        T = get_uniform_transition_matrix(prandom, nstates, 1)
        self.assertTrue(np.allclose(T, M))
        trans_a = UniformTransitionObject(prandom, nstates)
        trans_b = MatrixTransitionObject(M)
        stat_m = get_stationary_distribution(M)
        stat_a = trans_a.get_stationary_distribution()
        stat_b = trans_b.get_stationary_distribution()
        self.assertTrue(np.allclose(stat_m, stat_a))
        self.assertTrue(np.allclose(stat_m, stat_b))

    def test_uniform_transition_matrix(self):
        prandom = .1
        nstates = 4
        T = get_uniform_transition_matrix(prandom, nstates, 1)
        self.assertTrue(np.allclose(np.sum(T, 0), np.ones(nstates)))
        self.assertTrue(np.allclose(np.sum(T, 1), np.ones(nstates)))
        expected = np.linalg.matrix_power(T, 10)
        observed = get_uniform_transition_matrix(prandom, nstates, 10)
        self.assertTrue(np.allclose(expected, observed))

    def test_uniform_transition_object_transition(self):
        prandom = .1
        nstates = 4
        T = get_uniform_transition_matrix(prandom, nstates, 1)
        expected = np.linalg.matrix_power(T, 10)
        f = UniformTransitionObject(prandom, nstates)
        observed = [[f.get_transition_probability(i, j, 10) for i in range(nstates)] for j in range(nstates)]
        self.assertTrue(np.allclose(expected, observed))

    def test_uniform_transition_object_distribution(self):
        prandom = .1
        nstates = 4
        f = UniformTransitionObject(prandom, nstates)
        observed = [f.get_stationary_probability(i) for i in range(nstates)]
        expected = np.ones(nstates) / float(nstates)
        self.assertTrue(np.allclose(observed, expected))


if __name__ == '__main__':
    unittest.main()
