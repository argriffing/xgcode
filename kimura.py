"""
Diffusion approximations of continuous Wright-Fisher population genetics.

This should eventually include approximations for quantities
such as sojourn time densities,
expected time until absorption,
and fixation probabilities.
Initially the focus will be on genic selection
in a diploid population with random mating
and free recombination (sites evolve independently).
Throughout this documentation the preferred and unpreferred alleles
are named with respect to positive and negative s respectively.
The ratio of fitnesses of preferred to unpreferred alleles is 1 / (1-s) .
The formulas mostly follow
On the Probability of Fixation of Mutant Genes in a Population
by Motoo Kimura 1962, except in this module we use selection s as the
difference between homozygote selection instead of between allele selection.
We will also refer the the 1968 paper
THE AVERAGE NUMBER OF GENERATIONS UNTIL FIXATION OF
A MUTANT GENE IN A FINITE POPULATION
by Motoo Kimura and Tomoko Ohta.
"""

import unittest
import math

import numpy as np
from scipy import linalg

import MatrixUtil
import wfengine

def G(N_diploid, s):
    N = N_diploid * 2
    return math.exp(-N*s)

def get_pfix_approx(N_diploid, s):
    """
    This is an approximation of the fixation probability.
    It is derived by applying equation (3) in Kimura 1962,
    where p is the proportion of the preferred gene in the population,
    which in this case is 1/(2N) because it is a new mutant.
    It is given directly as equation (10).
    """
    N = N_diploid * 2
    if s:
        return math.expm1(-s) / math.expm1(-N*s)
    else:
        return 1.0 / N

class TestKimura(unittest.TestCase):
    def test_exact_fixation_probability(self):
        s = 0.03
        N_diploid = 3
        N_haploid = N_diploid * 2
        # Define a transition matrix with boundary conditions.
        # This approach is nice because the computationally intensive part
        # is just finding the equilibrium distribution of a process,
        # and you also get the conditional equilibrium distribution
        # of the dimorphic states.
        P = np.exp(wfengine.create_genic_diallelic(N_diploid, s))
        P[0, 0] = 0
        P[0, 1] = 1
        P[N_haploid, N_haploid] = 0
        P[N_haploid, 1] = 1
        v = MatrixUtil.get_stationary_distribution(P)
        fpa = v[-1] / (v[0] + v[-1])
        # Define a transition matrix with boundary conditions.
        # This approach is nice because it gives fixation probabilities
        # conditional on each possible initial allele frequency.
        P = np.exp(wfengine.create_genic_diallelic(N_diploid, s))
        A = P - np.eye(N_haploid + 1)
        b = np.zeros(N_haploid + 1)
        A[0, 0] = 1
        A[N_haploid, N_haploid] = 1
        b[0] = 0
        b[N_haploid] = 1
        x = linalg.solve(A, b)
        fpb = x[1]
        # Compare the transition probabilities.
        print fpa
        print get_pfix_approx(N_diploid, s)
        self.assertTrue(np.allclose(fpa, fpb))

if __name__ == '__main__':
    unittest.main()

