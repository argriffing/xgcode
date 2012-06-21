"""
Population genetics Markov models.

The model is similar to that of Wright-Fisher
and has mutation, selection, and recombination.
Each of these effects is associated with a separate global parameter.
For the initial and final state fields,
each row corresponds to a single chromosome
and each column corresponds to a position within the chromosome.
At every position in every chromosome
is either a low fitness or a high fitness allele,
denoted by 0 or 1 respectively.
"""

import unittest

import numpy as np
import gmpy

def bin_to_int(arr):
    """
    First entries of the input sequence are highest order.
    @param arr: sequence of zeros and ones
    @return: a nonnegative python integer
    """
    x = sum(v<<i for i, v in enumerate(reversed(arr)))
    return int(x)

def int_to_bin(x, npositions):
    """
    First entries of the output sequence are highest order.
    @param x: a nonnegative python integer
    @param npositions: length of the output array
    @return: an ndarray of zeros and ones
    """
    K = np.array(
            [gmpy.getbit(x, npositions-1-i) for i in range(npositions)],
            dtype=np.int8)
    return K

def bin2d_to_int(K):
    """
    @param K: an ndarray of integer ones and zeros
    @return: a python integer representing the state
    """
    return bin_to_int(np.ravel(K))

def int_to_bin2d(x, nrows, ncols):
    """
    @param x: a nonnegative python integer
    @param nrows: number of rows (chromosomes, domain-specifically)
    @param ncols: number of columns (positions, domain-specifically)
    @return: an ndarray with shape (nrows, ncols) containing zeros and ones
    """
    return int_to_bin(x, nrows*ncols).reshape(nrows, ncols)

def bitphase_to_nchanges(bitphase, npositions):
    """
    @param bitphase: a python integer
    @param npositions: length of binary array represented by the bitphase
    @return: the number of state changes along the binary array
    """
    nboundaries = npositions - 1
    return gmpy.popcount(gmpy.lowbits(bitphase ^ (bitphase >> 1), nboundaries))


def get_chromosome_distn(selection, recombination, K):
    """
    Define the distribution over child chromosomes.
    Given a parental population,
    the child chromosomes are independently and identically distributed
    when only selection and recombination are considered.
    This transition is independent from the immediately
    subsequent mutation effect.
    @param selection: a fitness ratio
    @param recombination: a linkage phase randomization probability
    @param K: ndarray parental population state
    @return: array defining a conditional distribution over chromosomes
    """
    nchromosomes, npositions = K.shape
    distn = np.zeros(1<<npositions)
    # sum over all ways to independently pick parental chromosomes
    # this is (nchromosomes)^2 things because repetition is allowed
    for parent_a in range(nchromosomes):
        weight_a = selection**np.sum(K[parent_a])
        for parent_b in range(nchromosomes):
            weight_b = selection**np.sum(K[parent_b])
            # sum over all recombination phases
            # this is 2^(npositions) things
            parent_pairs = zip(K[parent_a], K[parent_b])
            for bitphase in range(1<<npositions):
                # count the number of phase transitions in the bitphase
                nchanges = bitphase_to_nchanges(bitphase, npositions)
                # compute the weight corresponding to this phase transition
                weight_phase = 0
                weight_phase += recombination**nchanges
                weight_phase += (1-recombination)**(npositions-1-nchanges)
                # get the corresponding chromosome index
                phases = int_to_bin(bitphase, npositions)
                w = [pair[phase] for phase, pair in zip(phases, parent_pairs)]
                index = bin_to_int(w)
                distn[index] += weight_a * weight_b * weight_phase
    return distn / np.sum(distn)


class TestPopGenMarkov(unittest.TestCase):

    def test_bin_to_int_a(self):
        b = (0, 0, 1, 1, 0)
        expected = 6
        self.assertEqual(expected, bin_to_int(b))

    def test_bin_to_int_b(self):
        b = (1, 0, 0, 1, 0, 1)
        expected = 37
        self.assertEqual(expected, bin_to_int(b))

    def test_int_to_bin_a(self):
        x = 6
        npositions = 5
        expected = [0, 0, 1, 1, 0]
        self.assertEqual(expected, int_to_bin(x, npositions).tolist())

    def test_int_to_bin_b(self):
        x = 37
        npositions = 6
        expected = [1, 0, 0, 1, 0, 1]
        self.assertEqual(expected, int_to_bin(x, npositions).tolist())

    def test_bitphase_changes_a(self):
        bitphase = bin_to_int((1, 1, 0, 0, 0, 1))
        npositions = 6
        expected = 2
        observed = bitphase_to_nchanges(bitphase, npositions)
        self.assertEqual(expected, observed)

    def test_bitphase_changes_b(self):
        bitphase = bin_to_int((1, 1, 0, 0, 0, 1))
        npositions = 7
        expected = 3
        observed = bitphase_to_nchanges(bitphase, npositions)
        self.assertEqual(expected, observed)

    def test_bin2d_to_int(self):
        K = np.array([
            [0, 0, 1],
            [1, 0, 1]])
        expected = 13
        observed = bin2d_to_int(K)
        self.assertEqual(expected, observed)

    def test_int_to_bin2d(self):
        x = 13
        expected = np.array([
            [0, 0, 1],
            [1, 0, 1]])
        observed = int_to_bin2d(x, 2, 3)
        self.assertTrue(np.allclose(expected, observed))



if __name__ == '__main__':
    unittest.main()

