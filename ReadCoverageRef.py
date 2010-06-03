"""
This is supposed to model observations with reference sequences.
In this module, observations are 4-tuples.
The first element is the reference base read count.
The remaining three elements are counts of the three other non-ref bases.
"""

import math
import unittest

import numpy as np

import Util
import StatsUtil
import ReadCoverage
import ReadCoverageGap
import DGRP

class WrappedPoisson:
    def __init__(self, expectation):
        self.expectation = expectation
    def sample_observation(self):
        return np.random.poisson(self.expectation)
    def get_log_likelihood(self, obs):
        return StatsUtil.poisson_log_pmf(obs, self.expectation)
    def get_likelihood(self, obs):
        return math.exp(StatsUtil.poisson_log_pmf(obs, self.expectation))

class GoodMultiCoverage(ReadCoverage.UniformMixture):
    """
    This is a distribution over total coverages.
    """
    def __init__(self, nominal_coverage, kmulticoverages):
        states = [WrappedPoisson(nominal_coverage*(i+1))
                for i in range(kmulticoverages)]
        ReadCoverage.UniformMixture.__init__(self, states)

def gen_RR_distns(r):
    """
    Yield distributions over the elements of an observation.
    @param r: sequencing randomization rate
    """
    d = [r/4.0]*4
    d[0] = 1 - 3*r/4.0
    yield d

def gen_RA_distns(r):
    """
    Yield distributions over the elements of an observation.
    @param r: sequencing randomization rate
    """
    for i in range(1, 4):
        d = [r/4.0]*4
        d[0] = .5 - r/4.0
        d[i] = .5 - r/4.0
        yield d

def gen_AA_distns(r):
    """
    Yield distributions over the elements of an observation.
    @param r: sequencing randomization rate
    """
    for i in range(1, 4):
        d = [r/4.0]*4
        d[i] = 1 - 3*r/4.0
        yield d

def gen_AB_distns(r):
    """
    Yield distributions over the elements of an observation.
    @param r: sequencing randomization rate
    """
    for i in range(1, 4):
        for j in range(i+1, 4):
            d = [r/4.0]*4
            d[i] = .5 - r/4.0
            d[j] = .5 - r/4.0
            yield d

class SinglePatternState:
    def __init__(self, distribution, coverage_distribution):
        """
        @param distribution: expected nucleotide distribution
        @param coverage_distribution: a GoodMultiCoverage object
        """
        self.distribution = distribution
        self.coverage_distribution = coverage_distribution
    def sample_observation(self):
        n = self.coverage_distribution.sample_observation()
        return np.random.multinomial(n, self.distribution)
    def get_log_likelihood(self, obs):
        n = sum(obs)
        accum = 0
        accum += self.coverage_distribution.get_log_likelihood(n)
        accum += StatsUtil.multinomial_log_pmf(self.distribution, obs)
        return accum
    def get_likelihood(self, obs):
        return math.exp(self.get_log_likelihood(obs))

class GoodState(ReadCoverage.Mixture):
    def __init__(self, dref, dchild, seqerr, nomcoverage, kmulticoverages):
        """
        @param dref: a branch length parameter
        @param dchild: a branch length parameter
        @param seqerr: probability of sequence randomization
        @param nomcoverage: nominal coverage
        @param kmulticoverages: allowed multiples of nominal coverage
        """
        mcov = GoodMultiCoverage(nomcoverage, kmulticoverages)
        # define the states
        r = seqerr
        RR_states = [SinglePatternState(d, mcov) for d in gen_RR_distns(r)]
        RA_states = [SinglePatternState(d, mcov) for d in gen_RA_distns(r)]
        AA_states = [SinglePatternState(d, mcov) for d in gen_AA_distns(r)]
        AB_states = [SinglePatternState(d, mcov) for d in gen_AB_distns(r)]
        # define the distributions
        RR = ReadCoverage.UniformMixture(RR_states)
        RA = ReadCoverage.UniformMixture(RA_states)
        AA = ReadCoverage.UniformMixture(AA_states)
        AB = ReadCoverage.UniformMixture(AB_states)
        states = (RR, RA, AA, AB)
        zygo_distn = DGRP.get_zygosity_distribution(dref, dchild)
        ReadCoverage.Mixture.__init__(self, states, zygo_distn)

class HMMRecent(GoodState):
    """
    A predominantly homozygous region.
    """
    def __init__(self, x, y, z, seqerr, nomcoverage, kmulticoverages):
        GoodState.__init__(self, x+y, z, seqerr, nomcoverage, kmulticoverages)

class HMMAncient(GoodState):
    """
    This region has many heterozygous states.
    """
    def __init__(self, x, y, z, seqerr, nomcoverage, kmulticoverages):
        GoodState.__init__(self, x, y+z, seqerr, nomcoverage, kmulticoverages)

class HMMGarbage(ReadCoverage.UniformMixture):
    """
    This region has states with ill-defined zygosity.
    """
    def __init__(self, low, med, high):
        states = [ReadCoverageGap.FlatState(4, x) for x in (low, med, high)]
        ReadCoverage.UniformMixture.__init__(self, states)


class TestReadCoverageRef:

    def test_gen_xx_distns(self):
        for r in (.001, .01, .1, .9):
            RR = list(gen_RR_distns(r))
            RA = list(gen_RA_distns(r))
            AA = list(gen_AA_distns(r))
            AB = list(gen_AB_distns(r))
            self.assertEqual(len(RR), 1)
            self.assertEqual(len(RA), 3)
            self.assertEqual(len(AA), 3)
            self.assertEqual(len(AB), 6)
            all_distns = RR + RA + AA + AB
            for d in all_distns:
                self.assertAlmostEqual(sum(d), 1.0)


if __name__ == '__main__':
    unittest.main()
