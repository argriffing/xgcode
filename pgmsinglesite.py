"""
Single site population genetic markov model.

This is for large populations where the genome is only one nucleotide long.
"""

import unittest
import math

import numpy as np

import StatsUtil
import MatrixUtil

def create_drift_selection_transition_matrix(npop, selection_ratio):
    """
    The states are indexed by the number of mutants.
    @param npop: total population size
    @param selection_ratio: a value larger than unity means mutants are fitter
    @return: a transition matrix
    """
    nstates = npop + 1
    P = np.zeros((nstates, nstates))
    for a in range(nstates):
        # compute the i.i.d probability of picking a mutant
        p = (selection_ratio * a) / (selection_ratio * a + (npop-a))
        for b in range(nstates):
            # These are from a binomial distribution
            # with npop trials and p probability of success per trial.
            # (n choose k) p^k (1-p)^(n-k)
            observed_n = b
            max_n = npop
            p_success = p
            P[a, b] = math.exp(StatsUtil.binomial_log_pmf(
                observed_n, max_n, p_success))
    return P

def create_mutation_transition_matrix(npop, mutation_ab, mutation_ba):
    """
    The states are indexed by the number of mutants.
    @param npop: total population size
    @param mutation_ab: wild-type to mutant transition probability
    @param mutation_ba: mutant to wild-type transition probability
    @return: a transition matrix
    """
    StatsUtil.assert_probability(mutation_ab)
    StatsUtil.assert_probability(mutation_ba)
    nstates = npop + 1
    P = np.zeros((nstates, nstates))
    for a in range(nstates):
        for n_mut_to_wild in range(a+1):
            ba_observed_n = n_mut_to_wild
            ba_max_n = a
            ba_p_success = mutation_ba
            ba_log_p = StatsUtil.binomial_log_pmf(
                    ba_observed_n, ba_max_n, ba_p_success)
            for n_wild_to_mut in range(npop - a + 1):
                ab_observed_n = n_wild_to_mut
                ab_max_n = npop - a
                ab_p_success = mutation_ab
                ab_log_p = StatsUtil.binomial_log_pmf(
                        ab_observed_n, ab_max_n, ab_p_success)
                #
                p = math.exp(ba_log_p + ab_log_p)
                b = a + n_wild_to_mut - n_mut_to_wild
                P[a, b] += p
    return P

class Testpgmsinglesite(unittest.TestCase):
    def test_drift_selection(self):
        npop = 10
        selection_ratio = 2.0
        P = create_drift_selection_transition_matrix(
                npop, selection_ratio)
        MatrixUtil.assert_transition_matrix(P)
    def test_mutation(self):
        npop = 10
        mutation_ab = 0.1
        mutation_ba = 0.2
        P = create_mutation_transition_matrix(
                npop, mutation_ab, mutation_ba)
        MatrixUtil.assert_transition_matrix(P)

if __name__ == '__main__':
    unittest.main()
