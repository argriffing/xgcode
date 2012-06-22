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

import popgenmarkov
import MatrixUtil


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
        # After a pair of parental chromosomes are chosen
        # according to selection, one of them is passed at random to the
        # child with no recombination.
        child_distn = np.zeros(1<<npositions)
        for chra, chrb, p in popgenmarkov._gen_parental_triples(
                parent_chroms, selection):
            child_distn[chra] += 0.5 * p
            child_distn[chrb] += 0.5 * p
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

