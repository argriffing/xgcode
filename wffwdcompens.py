"""
Wright-Fisher forward simulation for compensatory mutations.

This model allows recombination and multiplicative genic selection.
The implementation in this module is slow,
and speeding it up the right way is tricky for unnecessary reasons.
Ideally it should be possible to speed it up by moving to cython,
but the problem is that cython does not have a nice interface for
random number sampling at C-like speeds; even in cython,
a call to a python/cython/numpy/scipy random number sampler will
destroy the speed of an inner loop by orders of magnitude.
The workaround that is available on the internet
is to avoid using the python/numpy/scipy random sampling functions altogether,
and instead to link to the gnu scientific library (gsl)
and use its random sampling functions instead.
.
This module is all about random sampling
so it is less straightforward to test.
.
As of the time of writing, np.random.choice only exists in
unreleased development versions of numpy (new in 1.7.0).
"""

import unittest
import random

import numpy as np

import wfcompens


AB_type = 0
Ab_type = 1
aB_type = 2
ab_type = 3

g_mutation = {
        (AB_type, 0) : aB_type,
        (AB_type, 1) : Ab_type,
        (Ab_type, 0) : ab_type,
        (Ab_type, 1) : AB_type,
        (aB_type, 0) : AB_type,
        (aB_type, 1) : ab_type,
        (ab_type, 0) : Ab_type,
        (ab_type, 1) : aB_type,
        }

def mutate_once(mutable_state, mutable_counts):
    """
    The two inputs are redundant for extra efficiency.
    The time complexity of this function should not depend on population size.
    @param mutable_state: sequence of haplotypes in the population
    @param mutable_counts: sequence of haplotype counts
    """
    # get the haploid population size
    N_hap = len(mutable_state)
    # pick the index of a random haploid in the population
    i = random.randrange(N_hap)
    # pick a site to mutate
    site = random.randrange(2)
    # Flip the state at the site.
    # This affects both the mutable_state and mutable_count arrays.
    # The haplotype coding is (AB, Ab, aB, ab).
    old_type = mutable_state[i]
    new_type = g_mutation[old_type, site]
    mutable_state[i] = new_type
    mutable_counts[old_type] -= 1
    mutable_counts[new_type] += 1
    # return None because this function works by mutating state
    return None

def recombine_once(mutable_state, mutable_counts):
    """
    The two inputs are redundant for extra efficiency.
    The time complexity of this function should not depend on population size.
    @param mutable_state: sequence of haplotypes in the population
    @param mutable_counts: sequence of haplotype counts
    """
    # get the haploid population size
    N_hap = len(mutable_state)
    # pick a random disjoint pair of haploid individuals from the population
    i, j = random.sample(xrange(N_hap), 2)
    # Change the population in-place according to recombination.
    # The haplotype coding is (AB, Ab, aB, ab).
    ti = mutable_state[i]
    tj = mutable_state[j]
    if (ti == AB_type and tj == ab_type) or (ti == ab_type and tj == AB_type):
        mutable_state[i] = Ab_type
        mutable_state[j] = aB_type
        mutable_counts[AB_type] -= 1
        mutable_counts[ab_type] -= 1
        mutable_counts[Ab_type] += 1
        mutable_counts[aB_type] += 1
    elif (ti == Ab_type and tj == aB_type) or (ti == aB_type and tj == Ab_type):
        mutable_state[i] = AB_type
        mutable_state[j] = ab_type
        mutable_counts[AB_type] += 1
        mutable_counts[ab_type] += 1
        mutable_counts[Ab_type] -= 1
        mutable_counts[aB_type] -= 1
    # return None because this function works by mutating state
    return None

def mutate_multiple(expected_events, mutable_state, mutable_counts):
    """
    @param expected_events: the expected number of mutation events
    @param mutable_state: sequence of haplotypes in the population
    @param mutable_counts: sequence of haplotype counts
    """
    nevents = np.random.poisson(expected_events)
    for i in range(nevents):
        mutate_once(mutable_state, mutable_counts)
    return None

def recombine_multiple(expected_events, mutable_state, mutable_counts):
    """
    @param expected_events: the expected number of mutation events
    @param mutable_state: sequence of haplotypes in the population
    @param mutable_counts: sequence of haplotype counts
    """
    nevents = np.random.poisson(expected_events)
    for i in range(nevents):
        recombine_once(mutable_state, mutable_counts)
    return None

def reselect(fitnesses, counts):
    """
    Get the allele counts for the new generation.
    Note that this uses haploid selection.
    When applied to diploid populations, this implies multiplicative selection.
    @param fitnesses: sequence of allele fitnesses
    @param counts: sequence of haploid allele counts
    @return: new counts
    """
    # get the haploid population size
    N_hap = np.sum(counts)
    # define the multinomial probabilities using multiplicative selection
    probs_kernel = np.array(fitnesses) * np.array(counts)
    probs = probs_kernel / np.sum(probs_kernel)
    # draw the new counts according to a multinomial distribution
    return np.random.multinomial(N_hap, probs)

def sample_hitting_time_slow(N_hap, mu, r, fitnesses):
    """
    This does not use the cython extension.
    Everything is with replacement.
    Sample a single hitting time in units of generational transitions.
    @param N_hap: haploid population size or twice diploid population size
    @param mu: expected number of mutations per generation
    @param r: expected number of recombinations per generation
    @param fitnesses: relative haploid fitnesses of haplotypes
    """
    ntransitions = 0
    counts = np.array([N_hap, 0, 0, 0])
    while not np.array_equal(counts, np.array([0, 0, 0, N_hap])):
        # the following line might be slow
        state = np.array(
                [0]*counts[0] + [1]*counts[1] + [2]*counts[2] + [3]*counts[3])
        mutate_multiple(mu, state, counts)
        recombine_multiple(r, state, counts)
        counts = reselect(fitnesses, counts)
        ntransitions += 1
    return ntransitions

def sample_hitting_time(N_hap, mu, r, fitnesses):
    """
    This uses the cython extension.
    This is currently with replacement,
    but should be modified to be without replacement
    when numpy 1.7.0 is released,
    or possibly when numpy 1.7.0 is packaged for ubuntu;
    this will probably in time for the 13.04 ubuntu release
    but not the 12.10 release.
    Sample a single hitting time in units of generational transitions.
    @param N_hap: haploid population size or twice diploid population size
    @param mu: expected number of mutations per generation
    @param r: expected number of recombinations per generation
    @param fitnesses: relative haploid fitnesses of haplotypes
    """
    ntransitions = 0
    haplotype_probs = np.zeros(4, dtype=float)
    counts = np.array([N_hap, 0, 0, 0], dtype=int)
    state = np.empty(N_hap, dtype=int)
    #target_counts = np.array([0, 0, 0, N_hap])
    #while not np.array_equal(counts, target_counts):
    while counts[-1] < N_hap:
        # expand the counts into the unary state
        wfcompens.expand_counts(state, counts)
        # do the mutation step
        nevents = np.random.poisson(mu)
        if nevents > N_hap:
            raise ValueError('sampling too many things without replacement')
        if nevents:
            #
            # TODO use np.random.choice when it becomes available
            #individuals = np.random.choice(
                    #N_hap, size=nevents, replace=False)
            individuals = np.random.randint(N_hap, size=nevents)
            #
            loci = np.random.randint(2, size=nevents)
            wfcompens.multiple_mutation(state, counts, individuals, loci)
        # do the recombination step
        nevents = np.random.poisson(r)
        if nevents*2 > N_hap:
            raise ValueError('sampling too many things without replacement')
        if nevents:
            #
            # TODO use np.random.choice when it becomes available
            #individuals = np.random.choice(
                    #N_hap, size=2*nevents, replace=False)
            individuals = np.random.randint(N_hap, size=2*nevents)
            #
            wfcompens.multiple_recombination(state, counts, individuals)
        # do the selection step
        wfcompens.reselection(haplotype_probs, fitnesses, counts)
        counts = np.random.multinomial(N_hap, haplotype_probs)
        # increment the number of observed generational transitions
        ntransitions += 1
    return ntransitions

def hitting_time_helper(N_hap, mu, r, s, nsamples):
    """
    @param N_hap: haploid population size
    @param mu: expected number of mutation events per generation
    @param r: expected number of recombination events per generation
    @param s: selection
    @param nsamples: average over this many samples
    """
    fitnesses = np.array([1, 1-s, 1-s, 1], dtype=float)
    total = 0
    for i in range(nsamples):
        g = sample_hitting_time(N_hap, mu, r, fitnesses)
        total += g
        print g
    return total / float(nsamples)

def test_hitting_time():
    """
    """
    N_diploid = 100
    N_hap = 2 * N_diploid
    theta = 0.01
    r = 0
    Ns = 0
    s = Ns / float(N_diploid)
    nsamples = 100
    mutation_rate = theta / 2
    generations = hitting_time_helper(N_hap, mutation_rate, r, s, nsamples)
    print 'theta:', theta
    print 'Ns:', Ns
    print 'generations:', generations
    print 'generations * theta:', generations * theta

def test_hitting_time_b():
    """
    I get 18,000 generations with mu = 0.8
    I get 8,400 generations with mu = 0.6
    I get 6,700 generations with mu = 0.5
    I get 5,800 generations with mu = 0.4
    I get 5,400 generations with mu = 0.3
    I get 5,800 generations with mu = 0.25
    I get 6,400 generations with mu = 0.2
    I get 8,700 generations with mu = 0.1
    I get 18,000 generations with mu = 0.05
    """
    N_diploid = 100
    N_hap = 2 * N_diploid
    mu = 0.25
    r = 0
    s = 0
    nsamples = 100
    ehit = hitting_time_helper(N_hap, mu, r, s, nsamples)
    print 'mu:', mu
    print 'expected hitting time:', ehit

def test_hitting_time_c():
    """
    This is about 7,000 generations mu = 0.5, Ns = 0.0
    This is about 4,600 generations mu = 0.5, Ns = 0.5
    This is about 4,500 generations mu = 0.5, Ns = 1.0
    This is about 5,300 generations mu = 0.5, Ns = 1.5
    This is about 7,200 generations mu = 0.5, Ns = 2.5
    """
    N_diploid = 100
    N_hap = 2 * N_diploid
    mu = 0.5
    r = 0
    Ns = 0
    s = Ns / float(N_diploid)
    nsamples = 100
    ehit = hitting_time_helper(N_hap, mu, r, s, nsamples)
    print 'mu:', mu
    print 'expected hitting time:', ehit

class TestForwardSampling(unittest.TestCase):
    def test_hitting_time(self):
        test_hitting_time()

if __name__ == '__main__':
    unittest.main()

