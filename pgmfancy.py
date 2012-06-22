"""
This is an extension of the popgenmarkov module.

Experimental features may include tradeoffs between
time and memory for matrix powers in endpoint conditioned path sampling,
unordering of chromosomes in a population,
and special functions for making inferences under population genetic
parameter edge cases, such as no mutation or no recombination.
The strange _s suffixes of functions in this module are to distinguish
functions that use the more complicated short state space
from the popgenmarkov functions that use the more straightforward
but less efficient longer state space.
"""

import unittest
from itertools import product
from collections import defaultdict

import numpy as np
import gmpy

import popgenmarkov
import MatrixUtil

def chroms_to_index(chroms, npositions):
    index = 0
    for chrom in chroms:
        index <<= npositions
        index |= chrom
    return index

def get_state_space_info(nchromosomes, npositions):
    """
    Get info related to the reduction from a large to a smaller state space.
    In this function ci means canonical index.
    @param nchromosomes: number of chromosomes in the population
    @param npositions: number of positions per chromosome
    @return: ci_to_short, short_to_count, sorted_chrom_lists
    """
    ci_to_short = {}
    short_to_count = []
    sorted_chrom_lists = []
    for chroms in product(range(1<<npositions), repeat=nchromosomes):
        sorted_chroms = sorted(chroms)
        # define the canonical long index
        ci = chroms_to_index(sorted_chroms, npositions)
        if ci not in ci_to_short:
            ci_to_short[ci] = len(sorted_chrom_lists)
            sorted_chrom_lists.append(sorted_chroms)
            short_to_count.append(1)
        else:
            short_to_count[ci_to_short[ci]] += 1
    return ci_to_short, short_to_count, sorted_chrom_lists

def get_mutation_transition_matrix_s(
        ci_to_short, short_to_count, sorted_chrom_lists,
        mutation, nchromosomes, npositions):
    nstates = len(sorted_chrom_lists)
    # map from ndiff to probability
    ndiff_to_p = np.zeros(nchromosomes * npositions + 1)
    for ndiff in range(nchromosomes * npositions + 1):
        nsame = nchromosomes * npositions - ndiff
        ndiff_to_p[ndiff] = (mutation**ndiff)*((1-mutation)**nsame)
    # define the mutation transition matrix
    P = np.zeros((nstates, nstates))
    for parent_chroms in sorted_chrom_lists:
        parent_ci = chroms_to_index(sorted(parent_chroms), npositions)
        parent_short = ci_to_short[parent_ci]
        for child_chroms in product(range(1<<npositions), repeat=nchromosomes):
            child_ci = chroms_to_index(sorted(child_chroms), npositions)
            child_short = ci_to_short[child_ci]
            child_index = chroms_to_index(child_chroms, npositions)
            diff = gmpy.hamdist(parent_ci, child_index)
            P[parent_short, child_short] += ndiff_to_p[diff]
    return P

def get_selection_recombination_transition_matrix_s(
        ci_to_short, short_to_count, sorted_chrom_lists,
        selection, recombination, nchromosomes, npositions):
    nstates = len(sorted_chrom_lists)
    # precompute conditional child chromosome distributions
    conditional_child_distns = popgenmarkov.precompute_conditional_child_distns(
            recombination, npositions)
    # init the unnormalized transition matrix
    P = np.zeros((nstates, nstates))
    for parent_short, parent_chroms in enumerate(sorted_chrom_lists):
        # define the distribution over parental chromosome ordered pairs
        parental_triples = list(popgenmarkov._gen_parental_triples(
            parent_chroms, selection))
        # define the distribution over child chromosomes
        child_distn = np.zeros(1<<npositions, dtype=float)
        for chra, chrb, p_parents in parental_triples:
            child_distn += p_parents * conditional_child_distns[chra, chrb]
        # choose child chromosomes independently
        for child_short, child_chroms in enumerate(sorted_chrom_lists):
            p = 1
            for chrom in child_chroms:
                p *= child_distn[chrom]
            P[parent_short, child_short] = short_to_count[child_short] * p
    return P

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
        source_index = chroms_to_index(parent_chroms, npositions)
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

    def test_invariant_selection_recombination_transition_s(self):
        selection = 1.1
        recombination = 0.01
        nchromosomes = 3
        npositions = 2
        ci_to_short, short_to_count, sorted_chrom_lists = get_state_space_info(
                nchromosomes, npositions)
        P = get_selection_recombination_transition_matrix_s(
                ci_to_short, short_to_count, sorted_chrom_lists,
                selection, recombination, nchromosomes, npositions)
        MatrixUtil.assert_transition_matrix(P)

    def test_invariant_mutation_transition_s(self):
        mutation = 0.01
        nchromosomes = 3
        npositions = 2
        ci_to_short, short_to_count, sorted_chrom_lists = get_state_space_info(
                nchromosomes, npositions)
        P = get_mutation_transition_matrix_s(
                ci_to_short, short_to_count, sorted_chrom_lists,
                mutation, nchromosomes, npositions)
        MatrixUtil.assert_transition_matrix(P)


if __name__ == '__main__':
    unittest.main()

