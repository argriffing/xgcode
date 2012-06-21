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
    return sum(v<<i for i, v in enumerate(reversed(arr)))

def int_to_bin(x, npositions):
    """
    First entries of the output sequence are highest order.
    @param x: a nonnegative python integer
    @param npositions: length of the output array
    @return: an ndarray of zeros and ones
    """
    return np.array([gmpy.getbit(x, npositions-1-i) for i in range(npositions)])

def multiline_state_to_ndarray(multiline_state):
    arr = []
    for line in stripped_lines(multiline_state.splitlines()):
        row = []
        for s in line:
            v = int(s)
            if v not in (0, 1):
                raise ValueError('invalid allele')
            row.append(v)
        arr.append(row)
    return np.array(arr)

def ndarray_to_multiline_state(K):
    return '\n'.join(''.join(str(x) for x in r) for r in K)

def ndarray_to_integer_state(K):
    """
    @param K: a rectangular array of integer ones and zeros
    @return: a python integer representing the state
    """
    nchromosomes, npositions = K.shape
    w = []
    for i in range(nchromosomes):
        for j in range(npositions):
            w.append(K[i,j])
    return sum(v<<i for i, v in enumerate(w))

def integer_state_to_ndarray(k, nchromosomes, npositions):
    """
    @param k: a python integer representing the state
    @return: a rectangular array of integer ones and zeros
    """
    arr = []
    for i in range(nchromosomes):
        row = []
        for j in range(npositions):
            v = k & 1
            k >>= 1
            row.append(v)
        arr.append(row)
    return np.array(arr)

def bitphase_to_nchanges(bitphase, npositions):
    mask = (1<<(npositions-1)) - 1
    return gmpy.popcount(mask & (bitphase ^ (bitphase >> 1)))



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
                w = [pair[(bitphase>>i) & 1] for i, pair in enumerate(
                    parent_pairs)]
                index = sum(v<<i for i, v in enumerate(w))
                distn[index] += weight_a * weight_b * weight_phase
    return distn / np.sum(distn)


class TestPopGenMarkov(unittest.TestCase):

    def test_bin_to_int(self):
        b = (0, 0, 1, 1, 0)
        expected = 6
        self.assertEqual(expected, bin_to_int(b))

    def test_int_to_bin(self):
        x = 6
        npositions = 5
        expected = [0, 0, 1, 1, 0]
        self.assertEqual(expected, int_to_bin(x, npositions).tolist())

    def test_bitphase_changes(self):
        #bitphase_to_nchanges(bitphase, npositions)
        pass

if __name__ == '__main__':
    unittest.main()
