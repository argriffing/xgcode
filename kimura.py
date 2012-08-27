"""
Diffusion approximations of continuous Wright-Fisher population genetics.

This should eventually include approximations for quantities
such as sojourn time densities,
expected time until absorption,
and fixation probabilities.
Initially the focus will be on genic selection
in a diploid population with random mating
and free recombination (sites evolve independently).
"""

import unittest

import numpy as np
from scipy import linalg

import MatrixUtil
import wfengine


class TestKimura(unittest.TestCase):
    def test_fixation_probability(self):
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
        self.assertTrue(np.allclose(fpa, fpb))

if __name__ == '__main__':
    unittest.main()

