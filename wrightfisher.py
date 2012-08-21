"""
Wright-Fisher transitions and their diffusion approximations.
"""

import unittest

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
    def test_foo(self):
        pass

if __name__ == '__main__':
    unittest.main()

