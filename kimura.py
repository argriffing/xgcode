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
The Kimura approximation of Fixation probability
with recessiveness/dominance h parameter can be done 
using wolfram alpha.
( integral from 0 to p of exp( a*x*(1 + h*(1-x)) ) ) ) /
( integral from 0 to 1 of exp( a*x*(1 + h*(1-x)) ) ) )
"""

import unittest
import math

import numpy as np
from scipy import linalg
from scipy import special

import MatrixUtil
import wfengine

def erfi(x):
    """
    Why does scipy.special not include this?
    """
    return 2 * x * special.hyp1f1(0.5, 1.5, x*x) / math.sqrt(math.pi)

def get_fixation_probability_chen(p, s, h):
    """
    This uses the parameter conventions from Christina Chen et al. 2008.
    @param p: initial allele frequency
    @param s: positive when the mutant allele is fitter
    @param h: dominance
    @return: fixation probability
    """
    s_eff = 2.0 * s
    beta = 2.0 * h - 1.0
    if not s_eff:
        return p
    if not beta:
        return math.expm1(-s_eff*p) / math.expm1(-s_eff)
    alpha = h / beta
    if beta * s_eff > 0:
        # overdominant if 0 < alpha < 1
        f = erfi
    elif beta * s_eff < 0:
        # underdominant if 0 < alpha < 1
        f = special.erf
    else:
        raise ValueError
    L = math.sqrt(abs(beta * s_eff))
    a0 = f(L*(0 - alpha))
    a1 = f(L*(p - alpha))
    b0 = f(L*(0 - alpha))
    b1 = f(L*(1 - alpha))
    pfix = (a1 - a0) / (b1 - b0)
    return pfix

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

def get_fixation_probabilities(P):
    """
    @param P: a wright fisher transition matrix
    """
    N_haploid = len(P) - 1
    A = P - np.eye(N_haploid + 1)
    b = np.zeros(N_haploid + 1)
    A[0, 0] = 1
    A[N_haploid, N_haploid] = 1
    b[0] = 0
    b[N_haploid] = 1
    return linalg.solve(A, b)

def get_fixation_conditioned_matrix(P, x):
    """
    UNFINISHED
    @param P: unconditioned transition matrix
    @param x: fixation probabilities
    @return: a smaller transition matrix
    """
    MatrixUtil.assert_transition_matrix(P)
    N_haploid = len(P) - 1
    B = A - np.eye(N_haploid)
    b = -np.ones(N_haploid)
    B[-1] = np.zeros(N_haploid)
    B[-1, -1] = 1
    b[-1] = 0
    y = linalg.solve(B, b)
    print y[0]

def J1_indefinite_integral(x, a):
    """
    Here is a wolfram alpha command that finds a limit of the indefinite
    integral of the integrand in J1 in equation (17) of Kimura and Ohto.
    limit x->1- of
    Ei(ax-a) - Ei(-ax) + log(x/(1-x)) -
    exp(-a)*( Ei(ax) - Ei(a-ax) + log((1-x)/x) )
    """
    eulergamma = -special.digamma(1)
    if x == 1:
        if a == 0:
            return 0
        else:
            x1 = math.exp(-a)
            x2 = math.log(abs(a)) + eulergamma
            x3 = math.exp(a) + 1
            x4 = special.expi(-a) + math.exp(-a) * special.expi(a)
            return x1 * x2 * x3 - x4
    else:
        x1 = special.expi(a*(x-1)) - special.expi(-a*x) + math.log(x/(1-x))
        x2 = special.expi(-a*(x-1)) - special.expi(a*x) + math.log(x/(1-x))
        return x1 + math.exp(-a)*x2

def J2_indefinite_integral(x, a):
    """
    limit x->0+ of
    Ei(ax) + Ei(-ax) - exp(a)*Ei(ax-a) - exp(-a)*Ei(a-ax) - 2*log(x/(1-x))
    """
    eulergamma = -special.digamma(1)
    if x == 0:
        if a == 0:
            # limit a->0 of
            # 2*(log|a| + eulergamma) - exp(a)*Ei(-a) - exp(-a)*Ei(a)
            return 0
        else:
            x1 = math.log(abs(a)) + eulergamma
            x2 = math.exp(a)*special.expi(-a)
            x3 = math.exp(-a)*special.expi(a)
            return 2*x1 - x2 - x3
    else:
        x1 = special.expi(a*x) + special.expi(-a*x)
        x2 = math.exp(a)*special.expi(a*(x-1))
        x3 = math.exp(-a)*special.expi(-a*(x-1))
        x4 = 2*math.log(x/(1-x))
        return x1 - x2 - x3 - x4

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
        #print fpa
        #print get_pfix_approx(N_diploid, s)
        self.assertTrue(np.allclose(fpa, fpb))

    def test_time_until_fixation_no_selection(self):
        s = 0
        N_diploid = 100
        N_haploid = N_diploid * 2
        print 'diploid population:', N_diploid
        P = np.exp(wfengine.create_genic_diallelic(N_diploid, s))
        # Create a transition matrix modified so that it is conditional
        # on eventual fixation.
        x = get_fixation_probabilities(P)
        print 'fixation probabilities:'
        print x
        A = (P * x)[1:, 1:]
        for i in range(N_haploid):
            A[i] /= np.sum(A[i])
        A[-1, -1] = 0
        A[-1, 0] = 1
        #
        #v = MatrixUtil.get_stationary_distribution(A)
        ## condition on a single mutation having already happened
        #v[0] -= v[-1]
        #v /= np.sum(v)
        #print 'expected generations until fixed, given eventual fixation:', (
                #(1 - v[-1]) / (v[-1]))
        #
        # Now use this conditional transition matrix
        # to set up a system of equations that will give
        # the expected number of generations until fixation
        # conditioned upon eventual fixation
        #
        B = A - np.eye(N_haploid)
        b = -np.ones(N_haploid)
        B[-1] = np.zeros(N_haploid)
        B[-1, -1] = 1
        b[-1] = 0
        y = linalg.solve(B, b)
        print 'expected time until fixation given some number of alleles:'
        print y
        #
        # Get the expected time until fixation given eventual fixation
        # when the population begins with proportion p mutants.
        #p = 5.0 * (1 / float(N_haploid))
        p = 0.1
        t0 = -2*N_haploid*(p/(1-p))*math.log(p)
        t1 = -(1/p)*2*N_haploid*(1-p)*math.log(1-p)
        print 'kimura initial p =', p
        print 'kimura expected to fixation:', t1
        print 'kimura expected to loss:', t0
        """
        # Create a transition matrix modified so that it is conditional
        # on eventual loss.
        A = (P * (1-x))[:-1, :-1]
        for i in range(N_haploid):
            A[i] /= np.sum(A[i])
        A[0, 0] = 0
        A[0, 1] = 1
        v = MatrixUtil.get_stationary_distribution(A)
        print 'expected generations until loss, given eventual loss:', (
                (1 - v[0]) / v[0])
        # Get the diffusion approximation
        p = 1 / float(N_haploid)
        t0 = -2*N_haploid*(p/(1-p))*math.log(p)
        print 'diffusion approximation:', t0
        """

    def test_expected_generations_until_absorption_no_selection(self):
        s = 0
        N_diploid = 100
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
        # Use the Kimura approximation to get the expected
        # number of generations until absorption for a neutral allele.
        p = 1 / float(N_haploid)
        #p_fixation = v[-1] / (v[0] + v[-1])
        #p_loss = v[0] / (v[0] + v[-1])
        # expected number of generations until fixation excluding loss
        t1 = -(1/p)*2*N_haploid*(1-p)*math.log(1-p)
        # expected number of generations until loss excluding fixation
        p_fixation = get_pfix_approx(N_diploid, s)
        p_loss = 1 - p_fixation
        t0 = -2*N_haploid*(p/(1-p))*math.log(p)
        t = p_loss * t0 + p_fixation * t1
        # foo
        #print t0
        #print t1
        #print t
        #print (1 - v[0] - v[-1]) / (v[0] + v[-1])

    def test_selection_coefficient_interpretations(self):
        for p in (0.1, 0.3, 0.6):
            for s in (-1.2, 0, 0.1, 1.5):
                # Compute the new genic probability
                # using the Zeng and Charlesworth model.
                f00 = 1.0
                f11 = 1.0 - s
                f01 = 0.5 * (f00 + f11)
                aa = f00 * p * p
                ab = f01 * p * (1-p)
                bb = f11 * (1-p) * (1-p)
                p_alpha = (aa + ab) / (aa + 2*ab + bb)
                # Compute using a difference
                # derived using wolfram alpha from the above model.
                delta = (0.5 * s) * p * (1 - p) / (1 + s*p - s)
                p_beta = p + delta
                # Check that the genic probabilities for the next generation
                # are the same using both formulas
                #print p_alpha
                #print p_beta
                self.assertTrue(np.allclose(p_alpha, p_beta))

    def test_fixation_time_expectation_symmetry(self):
        """
        The expected time to fixation given eventual fixation has a symmetry.
        If this is fix(p, s) then loss(p, s) = fix(1-p, -s).
        """
        pass

if __name__ == '__main__':
    unittest.main()

