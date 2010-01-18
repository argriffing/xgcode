"""
Clustering related algorithms.
"""

import StringIO
import unittest
import random

import numpy as np
import scipy

import Euclid
import NeighborJoining

# This distance matrix is from a neighbor joining paper by Lior Pachter.
g_lior = [
        [0.0, 3.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0],
        [3.0, 0.0, 3.0, 3.0, 3.0, 2.0, 2.0, 2.0],
        [2.0, 3.0, 0.0, 0.1, 0.4, 3.0, 3.0, 3.0],
        [2.0, 3.0, 0.1, 0.0, 0.4, 3.0, 3.0, 3.0],
        [2.0, 3.0, 0.4, 0.4, 0.0, 3.0, 3.0, 3.0],
        [3.0, 2.0, 3.0, 3.0, 3.0, 0.0, 0.1, 0.4],
        [3.0, 2.0, 3.0, 3.0, 3.0, 0.1, 0.0, 0.4],
        [3.0, 2.0, 3.0, 3.0, 3.0, 0.4, 0.4, 0.0]]

# This perturbed distance matrix is from a neighbor joining paper by Lior Pachter.
g_lior_perturbed = [
        [  0, 2.7, 2.6, 2.6, 2.6, 4.4, 4.4, 4.4],
        [2.7,   0, 4.4, 4.4, 4.4, 2.6, 2.6, 2.6],
        [2.6, 4.4,   0, 0.1, 0.4, 2.7, 2.7, 2.7],
        [2.6, 4.4, 0.1,   0, 0.4, 2.7, 2.7, 2.7],
        [2.6, 4.4, 0.4, 0.4,   0, 2.7, 2.7, 2.7],
        [4.4, 2.6, 2.7, 2.7, 2.7,   0, 0.1, 0.4],
        [4.4, 2.6, 2.7, 2.7, 2.7, 0.1,   0, 0.4],
        [4.4, 2.6, 2.7, 2.7, 2.7, 0.4, 0.4,   0]]

# This singular matrix is modified from a neighbor joining paper by Lior Pachter.
g_lior_singular = [
        [0.0, 3.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0],
        [3.0, 0.0, 3.0, 3.0, 3.0, 2.0, 2.0, 2.0],
        [2.0, 3.0, 0.0, 0.1, 0.4, 3.0, 3.0, 3.0],
        [2.0, 3.0, 0.1, 0.0, 0.4, 3.0, 3.0, 3.0],
        [2.0, 3.0, 0.4, 0.4, 0.0, 3.0, 3.0, 3.0],
        [3.0, 2.0, 3.0, 3.0, 3.0, 0.0, 0.1, 0.1],
        [3.0, 2.0, 3.0, 3.0, 3.0, 0.1, 0.0, 0.0],
        [3.0, 2.0, 3.0, 3.0, 3.0, 0.1, 0.0, 0.0]]

def get_R_stone(distance_matrix):
    """
    Return a function of the distance matrix (D^-1 - (D^-1)*1*1T*(D^-1)/(1T*D^-1*1)).
    The return value is negative one half of a graph Laplacian.
    @return: a function of the distance matrix
    """
    D = np.array(distance_matrix)
    n = len(D)
    D_inv = scipy.linalg.inv(D)
    ones = np.ones((n, n))
    numerator = np.dot(D_inv, np.dot(ones, D_inv))
    denominator = sum(sum(D_inv))
    R = D_inv - numerator / denominator
    return R

def get_R_stone_pinv(distance_matrix):
    """
    Return a function of the distance matrix (D^-1 - (D^-1)*1*1T*(D^-1)/(1T*D^-1*1)).
    Directly substitute the pseudoinverse for the inverse.
    The return value is negative one half of a graph Laplacian.
    @return: a function of the distance matrix
    """
    D = np.array(distance_matrix)
    n = len(D)
    D_pinv = scipy.linalg.pinv(D)
    ones = np.ones((n, n))
    numerator = np.dot(D_pinv, np.dot(ones, D_pinv))
    denominator = sum(sum(D_pinv))
    R = D_pinv - numerator / denominator
    return R

def get_R_balaji(distance_matrix):
    """
    Return a function of the distance matrix (D^-1 - (D^-1)*1*1T*(D^-1)/(1T*D^-1*1)).
    Use an equation I found in a paper by Balaji that accepts singular distance matrices.
    On Euclidean Distance Matrices.
    The return value is negative one half of a graph Laplacian.
    @return: a function of the distance matrix
    """
    D = np.array(distance_matrix)
    n = len(D)
    P = np.eye(n) - np.ones((n,n))/n 
    R_pinv = np.dot(P, np.dot(D, P)) 
    R = scipy.linalg.pinv(R_pinv)
    return R

def gen_assignments(n):
    """
    Generates Y vectors using the notation of Eric Stone.
    @param n: the number of states
    """
    if n < 2:
        raise ValueError('n must be at least two')
    arr = [1]*n
    while True:
        if set(arr) == set([-1, 1]):
            yield tuple(arr)
        for i in range(1, n):
            arr[i] *= -1
            if arr[i] == -1:
                break
        if arr == [1]*n:
            break

def get_exact_criterion(R_in, Y_in):
    """
    This slow function evaluates a partition of taxa.
    @param R_in: a particular transformation of the distance matrix
    @param Y_in: an array (row vector) with elements in {1, -1}
    @return: the value of the partition implied by Y
    """
    if not (set(Y_in) <= set([1, -1])):
        raise ValueError('Y_in must be a sequence of elements in {1, -1}')
    Y = np.array(Y_in)
    return np.dot(Y, np.dot(R_in, Y.T))


class DistanceMatrixSplitter:

    def get_complexity(self, n):
        """
        @param n: the number of rows of the square distance matrix
        @return: a value proportional to the expected computation time
        """
        raise NotImplementedError('override me')

    def _get_any_selection(self, distance_matrix):
        """
        @param distance_matrix: a numpy or row major distance matrix
        @return: a set of selected indices representing one of the two parts of the bipartition
        """
        raise NotImplementedError('override me')

    def get_selection(self, distance_matrix):
        """
        Using an algorithm specific to the derived class, find a bipartition.
        Of the two sets composing the partition, return the smaller.
        @param distance_matrix: a numpy or row major distance matrix
        @return: a set of indices
        """
        selection = self._get_any_selection(distance_matrix)
        n = len(distance_matrix)
        complement = set(range(n)) - selection
        if len(selection) < len(complement):
            return selection
        else:
            return complement


class UnnormalizedSpectralSignCutDMS(DistanceMatrixSplitter):

    def get_complexity(self, n):
        """
        @param n: the number of rows of the square distance matrix
        @return: a value proportional to the expected computation time
        """
        return n**3

    def _get_any_selection(self, D):
        """
        @param D: a numpy or row major distance matrix
        @return: a set of selected indices representing one of the two parts of the bipartition
        """
        # The fiedler eigenvector is calculated for the laplacian matrix associated with the distance matrix.
        # The signs of the elements of the fiedler eigenvector determine group assignment.
        L = Euclid.edm_to_laplacian(np.array(D))
        w, v = scipy.linalg.eigh(L)
        eigenvalue_info = list(sorted((abs(x), i) for i, x in enumerate(w)))
        stationary_eigenvector_index = eigenvalue_info[0][1]
        fiedler_eigenvector_index = eigenvalue_info[1][1]
        fiedler_eigenvector = v.T[fiedler_eigenvector_index]
        index_selection = set(i for i, value in enumerate(fiedler_eigenvector) if value > 0)
        return index_selection


class AffinityMinCutDMS(DistanceMatrixSplitter):
    """
    This method may fail if the laplacian matrix has a positive off-diagonal element.
    """

    def get_complexity(self, n):
        """
        @param n: the number of rows of the square distance matrix
        @return: a value proportional to the expected computation time
        """
        return n**3

    def _get_any_selection(self, D):
        """
        @param D: a numpy or row major distance matrix
        @return: a set of selected indices representing one of the two parts of the bipartition
        """
        # Get the selection corresponding to the min cut of the affinity matrix associated with the distance matrix.
        # Note that the negative off diagonals of the laplacian matrix form the affinity matrix.
        #L = Euclid.edm_to_laplacian(np.array(D))
        #index_selection = stoer_wagner_min_cut((-L).tolist())
        #return index_selection
        #TODO use brute force instead of the stoer wagner cut
        raise NotImplementedError()


class NeighborJoiningDMS(DistanceMatrixSplitter):

    def get_complexity(self, n):
        """
        @param n: the number of rows of the square distance matrix
        @return: a value proportional to the expected computation time
        """
        return n**2

    def _get_any_selection(self, distance_matrix):
        """
        @param distance_matrix: a numpy or row major distance matrix
        @return: a set of selected indices representing one of the two parts of the bipartition
        """
        return set(NeighborJoining.get_neighbors(distance_matrix))


class RandomDMS(DistanceMatrixSplitter):

    def get_complexity(self, n):
        """
        @param n: the number of rows of the square distance matrix
        @return: a value proportional to the expected computation time
        """
        return n

    def _get_any_selection(self, distance_matrix):
        """
        Make only informative splits, that is, the minimum selection size should be two.
        @param distance_matrix: a numpy or row major distance matrix
        @return: a set of selected indices representing one of the two parts of the bipartition
        """
        n = len(distance_matrix)
        min_selection_count = 2
        max_selection_count = n - 2
        if max_selection_count < min_selection_count:
            raise ValueError('the distance matrix is too small to make an informative split')
        selection_count = random.randrange(min_selection_count, max_selection_count + 1)
        index_selection = set(random.sample(range(n), selection_count))
        return index_selection

class StoneExactDMS(DistanceMatrixSplitter):

    def get_complexity(self, n):
        """
        @param n: the number of rows of the square distance matrix
        @return: a value proportional to the expected computation time
        """
        return (n**2) * (2**n)

    def _get_any_selection(self, distance_matrix):
        """
        @param distance_matrix: a numpy or row major distance matrix
        @return: a set of selected indices representing one of the two parts of the bipartition
        """
        n = len(distance_matrix)
        R = get_R_balaji(distance_matrix)
        # find the best Y vector, where each element of Y is one or negative one
        best_value_vector_pair = None
        for assignment in gen_assignments(n):
            Y = np.array(assignment)
            value = get_exact_criterion(R, Y)
            value_vector_pair = (value, Y)
            if (best_value_vector_pair is None) or (value_vector_pair[0] > best_value_vector_pair[0]):
                best_value_vector_pair = value_vector_pair
        # get the index set associated with the best vector
        best_value, best_vector = best_value_vector_pair
        index_selection = set(i for i, element in enumerate(best_vector) if element > 0)
        return index_selection


class StoneSpectralSignDMS(DistanceMatrixSplitter):

    def get_complexity(self, n):
        """
        @param n: the number of rows of the square distance matrix
        @return: a value proportional to the expected computation time
        """
        return n**3

    def _get_any_selection(self, distance_matrix):
        """
        @param distance_matrix: a numpy or row major distance matrix
        @return: a set of selected indices representing one of the two parts of the bipartition
        """
        n = len(distance_matrix)
        R = get_R_balaji(distance_matrix)
        # The signs of the elements of the fiedler eigenvector determine group assignment.
        w, v = scipy.linalg.eigh(R)
        eigenvalue_info = list(sorted((abs(x), i) for i, x in enumerate(w)))
        stationary_eigenvector_index = eigenvalue_info[0][1]
        fiedler_eigenvector_index = eigenvalue_info[1][1]
        fiedler_eigenvector = v.T[fiedler_eigenvector_index]
        index_selection = set(i for i, value in enumerate(fiedler_eigenvector) if value > 0)
        return index_selection

class StoneSpectralThresholdDMS(DistanceMatrixSplitter):

    def get_complexity(self, n):
        """
        @param n: the number of rows of the square distance matrix
        @return: a value proportional to the expected computation time
        """
        return n**3

    def _get_any_selection(self, distance_matrix):
        """
        @param distance_matrix: a numpy or row major distance matrix
        @return: a set of selected indices representing one of the two parts of the bipartition
        """
        n = len(distance_matrix)
        R = get_R_balaji(distance_matrix)
        # define the Fielder eigenvector
        w, v = scipy.linalg.eigh(R)
        eigenvalue_info = list(sorted((abs(x), i) for i, x in enumerate(w)))
        stationary_eigenvector_index = eigenvalue_info[0][1]
        fiedler_eigenvector_index = eigenvalue_info[1][1]
        fiedler_eigenvector = v.T[fiedler_eigenvector_index]
        # find the bipartition defined by the best cut of the fiedler eigenvector
        element_index_pairs = [(element, i) for i, element in enumerate(fiedler_eigenvector)]
        sorted_indices = [i for element, i in sorted(element_index_pairs)]
        best_value_vector_pair = None
        for cut in range(n-1):
            assignment = [0]*n
            for i, index in enumerate(sorted_indices):
                if i < cut + 1:
                    assignment[index] = -1
                else:
                    assignment[index] = 1
            v = np.array(assignment)
            value = np.dot(v, np.dot(R, v.T))
            value_vector_pair = (value, v)
            if (best_value_vector_pair is None) or (value_vector_pair[0] > best_value_vector_pair[0]):
                best_value_vector_pair = value_vector_pair
        # get the index set associated with the best vector
        best_value, best_vector = best_value_vector_pair
        index_selection = set(i for i, element in enumerate(best_vector) if element < 0)
        return index_selection



class TestClustering(unittest.TestCase):

    def test_pseudoinverse_a(self):
        """
        A test comparing the Stone and Balaji equations.
        """
        D = g_lior
        self.assertTrue(np.allclose(get_R_stone(D), get_R_balaji(D)))

    def test_pseudoinverse_b(self):
        """
        A test comparing the Stone and Balaji equations.
        """
        D = g_lior_perturbed
        self.assertTrue(np.allclose(get_R_stone(D), get_R_balaji(D)))

    def test_pseudoinverse_c(self):
        """
        A test comparing the Stone and Balaji equations.
        """
        D = g_lior_singular
        try:
            stone = get_R_stone(D)
        except scipy.linalg.LinAlgError, e:
            pass
        else:
            self.fail('inverting a singular matrix is supposed to fail')
        try:
            balaji = get_R_balaji(D)
        except scipy.linalg.LinAlgError, e:
            self.fail('the pseudoinverse is supposed to work even on a singular matrix')

    def test_pseudoinverse_d(self):
        """
        A test comparing the Stone and Balaji equations.
        """
        D = g_lior_singular
        self.assertTrue(np.allclose(get_R_stone_pinv(D), get_R_balaji(D)))

    def test_exact_criterion(self):
        D = g_lior
        n = len(D)
        R = get_R_balaji(D)
        value_vector_pairs = [(get_exact_criterion(R, Y), Y) for Y in gen_assignments(n)]
        expected_count = (1<<(n-1)) - 1
        observed_count = len(value_vector_pairs)
        self.assertEqual(expected_count, observed_count)
        # the best spectral approximation should be no better than the best exact partition
        splitter = StoneSpectralSignDMS()
        spectral_index_set = splitter.get_selection(D)
        spectral_Y = np.array([(1 if i in spectral_index_set else -1) for i in range(n)])
        spectral_value = get_exact_criterion(R, spectral_Y)
        best_exact_value, best_exact_Y = max(value_vector_pairs)
        if spectral_value > best_exact_value:
            self.fail('the spectral approximation should be no better than the exact criterion')


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestClustering)
    unittest.TextTestRunner(verbosity=2).run(suite)
