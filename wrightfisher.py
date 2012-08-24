"""
Wright-Fisher transitions and their diffusion approximations.

This also includes testing for the wfengine cython module.
"""

import unittest

import numpy as np

import MatrixUtil
import wfengine


def genic_diallelic(fi, fj, ni, nj):
    """
    The two returned probabilities sum to 1.
    This assumes that the fitness of the heterozygote is the mean
    of the fitness of the homozygotes,
    and that the parent population has exactly hardy weinberg proportions.
    These are the assumptions made, for example, in the paper by Kai Zeng 2010.
    Note that ni + nj = 2N where N is the diploid population size.
    @param fi: fitness of allele i
    @param fj: fitness of allele j
    @param ni: count of allele i in parent population
    @param nj: count of allele j in parent population
    @return: (ci, cj) iid allele probabilities for the child population
    """
    # get the parent allele proportions
    pi = ni / float(ni + nj)
    pj = nj / float(ni + nj)
    # get the fitnesses of the heterozygote and both homozygotes
    fij = 0.5 * (fi + fj)
    fii = fi
    fjj = fj
    # get the child probabilities
    a = fii*pi*pi + fij*pi*pj
    b = fjj*pj*pj + fij*pi*pj
    ci = a / float(a+b)
    cj = b / float(a+b)
    return ci, cj
    

class TestWrightFisher(unittest.TestCase):
    def test_wright_fisher(self):
        N_diallelic = 3
        s = 0.03
        # compute the wright fisher transition matrix directly
        Ma = np.exp(wfengine.create_genic_diallelic(N_diallelic, s))
        MatrixUtil.assert_transition_matrix(Ma)
        # compute the wright fisher transition matrix less directly
        N = 2*N_diallelic
        fi = 1.0
        fj = 1.0 - s
        log_distns = np.zeros((N+1, 2))
        for h in range(0, N+1):
            probs = genic_diallelic(fi, fj, h, (N-h))
            log_distns[h] = np.log(np.array(probs))
        Mb = np.exp(wfengine.expand_multinomials(N, log_distns))
        MatrixUtil.assert_transition_matrix(Mb)
        # compare the transition matrices
        self.assertTrue(np.allclose(Ma, Mb))


if __name__ == '__main__':
    unittest.main()

