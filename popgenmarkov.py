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
from itertools import product

import numpy as np
import gmpy

import MatrixUtil

def bin_to_int(arr):
    """
    First entries of the input sequence are highest order.
    Speeding up this function would help a lot.
    @param arr: sequence of zeros and ones
    @return: a nonnegative python integer
    """
    """
    x = 0
    for v in arr:
        x = (x << 1) | v
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

#TODO this has been copied into another function and may be unused
def _parent_pair_to_child_distn(
        chrom_a, chrom_b, recombination, npositions, phase_distn):
    """
    In this function chromosomes are python integers.
    The alleles are defined by bits of the binary chromosome representation.
    The output of this function should be cached.
    @param chrom_a: first parent chromosome as a python integer
    @param chrom_b: second parent chromosome as a python integer
    @param recombination: a linkage phase change probability
    @param npositions: number of positions per chromosome
    @param phase_distn: a probability distribution over phases
    @return: a distribution over child chromosomes
    """
    # define the size of the chromosome state space
    n = 1 << npositions
    distn = np.zeros(n, dtype=float)
    for phase in range(n):
        #nchanges = bitphase_to_nchanges(phase, npositions)
        #p_phase = 0
        #p_phase += recombination**nchanges
        #p_phase += (1-recombination)**(npositions-1-nchanges)
        p_phase = phase_distn[phase]
        chrom_c = (chrom_a & phase) | (chrom_b & ~phase)
        distn[chrom_c] += p_phase
    return distn

def precompute_conditional_child_distns(recombination, npositions):
    """
    @param recombination: a linkage phase change probability
    @param npositions: number of positions per chromosome
    @return: a shape (1<<npositions, 1<<npositions, 1<<npositions) ndarray
    """
    # define the size of the chromosome state space
    n = 1 << npositions
    # precompute phase distribution
    phase_distn = np.zeros(n, dtype=float)
    for phase in range(n):
        nchanges = bitphase_to_nchanges(phase, npositions)
        p = 0.5
        p *= recombination**nchanges
        p *= (1-recombination)**(npositions-1-nchanges)
        phase_distn[phase] = p
    # precompute the conditional child distributions
    conditional_child_distns = np.zeros((n, n, n), dtype=float)
    for chra in range(n):
        for chrb in range(n):
            for phase in range(n):
                chrc = (chra & phase) | (chrb & ~phase)
                conditional_child_distns[chra, chrb, chrc] += phase_distn[phase]
    return conditional_child_distns

def _gen_parental_triples(population, selection):
    """
    This gives a sparse distribution over pairs of parental chromosomes.
    Yield a sequence of (chra, chrb, probability) triples.
    @param population: a sequence of chromosomes each as python integers
    @param selection: a fitness ratio
    """
    n = len(population)
    # get the distribution over indices into the parental population
    distn = np.zeros(n)
    for i, chrom in enumerate(population):
        distn[i] = selection**gmpy.popcount(chrom)
    distn /= np.sum(distn)
    # Define the triples assuming that parental chromosomes
    # are drawn independently according to relative fitness.
    for chra, pa in zip(population, distn):
        for chrb, pb in zip(population, distn):
            yield chra, chrb, pa*pb

def get_selection_recombination_transition_matrix(
        selection, recombination, nchromosomes, npositions):
    """
    """
    nstates = 1 << (nchromosomes * npositions)
    # precompute conditional child chromosome distributions
    conditional_child_distns = precompute_conditional_child_distns(
            recombination, npositions)
    # init the unnormalized transition matrix
    P = np.zeros((nstates, nstates))
    for parent_chroms in product(range(1<<npositions), repeat=nchromosomes):
        # define the source index
        source_index = 0
        for chrom in parent_chroms:
            source_index <<= npositions
            source_index |= chrom
        # define the distribution over parental chromosome ordered pairs
        parental_triples = list(_gen_parental_triples(parent_chroms, selection))
        # define the distribution over child chromosomes
        child_distn = np.zeros(1<<npositions, dtype=float)
        for chra, chrb, p_parents in parental_triples:
            child_distn += p_parents * conditional_child_distns[chra, chrb]
        # choose child chromosomes independently
        for child_chroms in product(range(1<<npositions), repeat=nchromosomes):
            # define the sink index and conditional probability
            p = 1
            sink_index = 0
            for chrom in child_chroms:
                p *= child_distn[chrom]
                sink_index <<= npositions
                sink_index |= chrom
            P[source_index, sink_index] = p
    return P

def get_chromosome_distn(selection, recombination, K):
    """
    Define the distribution over child chromosomes.
    Given a parental population,
    the child chromosomes are independently and identically distributed
    when only selection and recombination are considered.
    This transition is independent from the immediately
    subsequent mutation effect.
    @param selection: a fitness ratio
    @param recombination: a linkage phase change probability
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
                weight_phase = 1
                weight_phase *= recombination**nchanges
                weight_phase *= (1-recombination)**(npositions-1-nchanges)
                # get the corresponding chromosome index
                phases = int_to_bin(bitphase, npositions)
                w = [pair[phase] for phase, pair in zip(phases, parent_pairs)]
                index = bin_to_int(w)
                distn[index] += weight_a * weight_b * weight_phase
    return distn / np.sum(distn)

def get_chromosome_distn_fast(selection, recombination, K):
    """
    This is a faster version with more bitwise cleverness.
    """
    nchromosomes, npositions = K.shape
    chromos = [bin_to_int(row) for row in K]
    distn = np.zeros(1<<npositions)
    # sum over all ways to independently pick parental chromosomes
    # this is (nchromosomes)^2 things because repetition is allowed
    for chromo_a in chromos:
        weight_a = selection**gmpy.popcount(chromo_a)
        for chromo_b in chromos:
            weight_b = selection**gmpy.popcount(chromo_b)
            for phase in range(1<<npositions):
                nchanges = bitphase_to_nchanges(phase, npositions)
                weight_phase = 1
                weight_phase *= recombination**nchanges
                weight_phase *= (1-recombination)**(npositions-1-nchanges)
                chromo_c = (chromo_a & phase) | (chromo_b & ~phase)
                distn[chromo_c] += weight_a * weight_b * weight_phase
    return distn / np.sum(distn)

def get_mutation_transition_matrix(mutation, nchromosomes, npositions):
    nstates = 1 << (nchromosomes * npositions)
    # map from ndiff to probability
    ndiff_to_p = np.zeros(nchromosomes * npositions + 1)
    for ndiff in range(nchromosomes * npositions + 1):
        nsame = nchromosomes * npositions - ndiff
        ndiff_to_p[ndiff] = (mutation**ndiff)*((1-mutation)**nsame)
    # define the mutation transition matrix
    P = np.zeros((nstates, nstates))
    for source in range(nstates):
        for sink in range(nstates):
            P[source, sink] = ndiff_to_p[gmpy.hamdist(source, sink)]
    return P


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

    def test_regress_chromo_distn(self):
        selection = 2.0
        recombination = 0.01
        K = np.array([
            [0, 0, 1],
            [1, 0, 1]])
        da = get_chromosome_distn(selection, recombination, K)
        db = get_chromosome_distn_fast(selection, recombination, K)
        self.assertTrue(np.allclose(da, db))

    def test_invariant_selection_recombination_transition(self):
        selection = 1.1
        recombination = 0.01
        nchromosomes = 3
        npositions = 2
        P = get_selection_recombination_transition_matrix(
                selection, recombination, nchromosomes, npositions)
        MatrixUtil.assert_transition_matrix(P)

    def test_invariant_mutation_transition(self):
        mutation = 0.01
        nchromosomes = 3
        npositions = 2
        P = get_mutation_transition_matrix(mutation, nchromosomes, npositions)
        MatrixUtil.assert_transition_matrix(P)

if __name__ == '__main__':
    unittest.main()

