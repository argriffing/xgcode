"""
This is a helper module for dealing with multinomial states.
In this module I am defining a multinomial state
as a sequence of k nonnegative integers whose sum is N.
The state space is parameterized by N and k.
This module is expected in the API for the wfengine.pyx python module.
"""

import unittest

import numpy as np

import Util


def get_nstates(N, k):
    """
    This number is on the order of O( N^(k-1) ).
    """
    return Util.choose(N+k-1, k-1)

def gen_states(N, k):
    """
    Yield multinomial states.
    Each state is a list of length k and with sum N.
    The state ordering is one which simple to generate recursively.
    @param N: population size
    @param k: number of bins
    """
    if k == 1:
        yield [N]
    elif k == 2:
        for i in range(N+1):
            yield [i, N-i]
    else:
        for i in range(N+1):
            for suffix in gen_states(N-i, k-1):
                yield [i] + suffix

def gen_reversed_states(N, k):
    """
    Yield multinomial states.
    Each state is a list of length k and with sum N.
    The state ordering is one which simple to generate recursively.
    @param N: population size
    @param k: number of bins
    """
    if k == 1:
        yield [N]
    elif k == 2:
        for i in range(N+1):
            yield [N-i, i]
    else:
        for i in range(N+1):
            for prefix in gen_reversed_states(N-i, k-1):
                yield prefix + [i]

def get_sorted_states(N, k):
    """
    Population states are sorted by decreasing max entry count.
    This implies that the k states that have
    one entry equal to N and the rest of the states equal to 0
    will be yielded first.
    @param N: size of haploid population, or twice a diploid population size
    @param k: number of alleles, for example 2 or 4
    @return: M where M[i,j] is count of allele j for pop state index i
    """
    return np.array(sorted(
        gen_reversed_states(N, k), reverse=True, key=max))

def get_inverse_map(M):
    """
    The input M[i,j] is count of allele j for pop state index i.
    The output T[(i,j,...)] maps allele count tuple to pop state index
    @param M: multinomial state map
    @return: T
    """
    nstates, k = M.shape
    N = np.sum(M[0])
    inverse_map_shape = tuple([N+1]*k)
    T = -1 * np.ones(inverse_map_shape, dtype=int)
    for i, state in enumerate(M):
        T[tuple(state)] = i
    return T

class TestMultinomialState(unittest.TestCase):
    def test_numpy_indexing(self):
        cube = np.array([[[1, 2], [3, 4]], [[5, 6], [7, 8]]])
        self.assertEqual(cube[0,0,1], 2)
        self.assertEqual(cube[(0,0,1)], 2)
        self.assertTrue(np.array_equal(
            cube[np.array([0,0,1])],
            [[[1,2],[3,4]], [[1,2],[3,4]], [[5,6],[7,8]]],
            ))
    def test_state_count(self):
        N = 3
        k = 2
        self.assertEqual(len(list(gen_states(N,k))), get_nstates(N,k))
    def test_state_generation(self):
        N=3
        k=2
        expected = [ [0, 3], [1, 2], [2, 1], [3, 0] ]
        observed = list(gen_states(N, k))
        self.assertEqual(expected, observed)
    def test_sorted_state_generation(self):
        N=3
        k=2
        expected = [ [3, 0], [0, 3], [2, 1], [1, 2] ]
        observed = get_sorted_states(N, k)
        self.assertEqual(expected, observed.tolist())
    def test_inverse_state_map(self):
        N=3
        k=2
        M = get_sorted_states(N, k)
        M_expected = [ [3, 0], [0, 3], [2, 1], [1, 2] ]
        self.assertEqual(M_expected, M.tolist())
        T = get_inverse_map(np.array(M))
        self.assertEqual(T[3, 0], 0)
        self.assertEqual(T[0, 3], 1)
        self.assertEqual(T[2, 1], 2)
        self.assertEqual(T[1, 2], 3)

if __name__ == '__main__':
    unittest.main()

