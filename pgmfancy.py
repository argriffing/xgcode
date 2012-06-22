"""
This is an extension of the popgenmarkov module.

Experimental features may include tradeoffs between
time and memory for matrix powers in endpoint conditioned path sampling,
unordering of chromosomes in a population,
and special functions for making inferences under population genetic
parameter edge cases, such as no mutation or no recombination.
"""

import unittest
from itertools import product

import numpy as np
import gmpy

import popgenmarkov
import MatrixUtil

def get_state_space_info(nchromosomes, npositions):
    """
    """
    pass

def get_selection_transition_matrix(selection, nchromosomes, npositions):
    """
    Note that this includes only selection and not recombination or mutation.
    Therefore the transition matrix will be very sparse.
    @param selection: a fitness ratio
    @param nchromosomes: number of chromosomes in the population
    @param npositions: number of positions per chromosome
    """
    # init the unnormalized transition matrix
    nstates = 1 << (nchromosomes * npositions)
    P = np.zeros((nstates, nstates))
    for parent_chroms in product(range(1<<npositions), repeat=nchromosomes):
        # define the source index
        source_index = 0
        for chrom in parent_chroms:
            source_index <<= npositions
            source_index |= chrom
        # get the distribution over indices into the parental population
        parent_index_distn = np.zeros(nchromosomes)
        for i, chrom in enumerate(parent_chroms):
            parent_index_distn[i] = selection**gmpy.popcount(chrom)
        parent_index_distn /= np.sum(parent_index_distn)
        # choose child chromosomes independently
        for parent_idxs in product(range(nchromosomes), repeat=nchromosomes):
            # define the sink index and conditional probability
            p = 1
            sink_index = 0
            for i in parent_idxs:
                p *= parent_index_distn[i]
                child_chrom = parent_chroms[i]
                sink_index <<= npositions
                sink_index |= child_chrom
            P[source_index, sink_index] += p
    return P


class TestPGMFancy(unittest.TestCase):

    def test_invariant_selection_transition(self):
        selection = 1.1
        nchromosomes = 3
        npositions = 2
        P = get_selection_transition_matrix(
                selection, nchromosomes, npositions)
        MatrixUtil.assert_transition_matrix(P)

    def test_regress_selection_transition(self):
        selection = 1.1
        recombination = 0.0
        nchromosomes = 3
        npositions = 2
        Pa = get_selection_transition_matrix(
                selection, nchromosomes, npositions)
        Pb = popgenmarkov.get_selection_recombination_transition_matrix(
                selection, recombination, nchromosomes, npositions)
        self.assertTrue(np.allclose(Pa, Pb))



if __name__ == '__main__':
    unittest.main()

