"""
Check a small recessive unpreferred nucleotide model.

This toy model of molecular evolution will have two free parameters
plus empirical equilibrium distribution parameters.
The first free parameter scales the rate matrix --
this parameter is constrained to be positive,
but this constraint is enforced by a log transformation.
The second free parameter determines the degree to which the unpreferred
allele is recessive --
this parameter is unconstrained,
and includes 'genic' selection, pure recessivity, pure dominance,
over-recessivity, over-dominance, and all other recessivity/dominance values.
"""

from StringIO import StringIO
import math
import cmath
import random

import numpy
import scipy
import scipy.special
import algopy

import Form
import FormOut

def kimura_1957_54_numerator_analytic(p, c, d):
    """
    From Kimura 1957, equation (5.4).
    This is the diffusion approximation of the fixation probability
    of an allele that has reached frequency p in the population,
    with scaled selection c = Ns
    and dominance/recessivity parameter d = 2h - 1.
    @param p: initial allele frequency in the population
    @param c: population-scaled selection coefficient
    @param d: transformed dominance/recessivity parameter
    @return: diffusion estimate of fixation probability
    """
    # Mathematica notation:
    # Integrate[Exp[-2*c*d*x*(1-x) - 2*c*x], {x, 0, p}]
    #
    # precompute some intermediate quantities
    sqrt_c = algopy.sqrt(c)
    sqrt_d = algopy.sqrt(d)
    sqrt_2 = math.sqrt(2)
    sqrt_pi = math.sqrt(math.pi)
    #
    # compute the numerator
    erfi_num_a = algopy.special.erfi(
            (sqrt_c * (1 + d)) / (sqrt_2 * sqrt_d))
    erfi_num_b = algopy.special.erfi(
            (sqrt_c * (1 + d - 2*d*p)) / (sqrt_2 * sqrt_d))
    num = (sqrt_pi / sqrt_2) * (erfi_num_a + erfi_num_b)
    #
    # compute the denominator
    exp_den_a = algopy.exp((c*((1+d)**2)) / (2*d))
    den = 2*algopy.sqrt(c)*algopy.sqrt(d)*exp_den_a
    #
    return num / den

def kimura_1957_54_denominator_analytic(c, d):
    """
    See the corresponding description for the numerator.
    This function computes the normalizing constant.
    """
    # Mathematica notation:
    # Integrate[Exp[-2*c*d*x*(1-x) - 2*c*x], {x, 0, 1}]
    #
    # precompute some intermediate quantities
    # FIXME: algopy has no csqrt but then again also does not support complex.
    #sqrt_c = cmath.sqrt(c)
    #sqrt_d = cmath.sqrt(d)
    sqrt_c = algopy.sqrt(c)
    sqrt_d = algopy.sqrt(d)
    sqrt_2 = math.sqrt(2)
    exp_2c = algopy.exp(2*c)
    #
    # compute the numerator
    dawsn_num_a = algopy.special.dawsn(
            (sqrt_c * (d-1)) / (sqrt_2 * sqrt_d))
    dawsn_num_b = algopy.special.dawsn(
            (sqrt_c * (d+1)) / (sqrt_2 * sqrt_d))
    num = dawsn_num_a + exp_2c * dawsn_num_b
    #
    # compute the denominator
    den = sqrt_2 * sqrt_c * sqrt_d * exp_2c
    #
    return num / den

def kimura_1957_54_denominator_analytic_b(c, d):
    """
    The numerator of (5.4) goes to 1 in the scaling limit.
    So we only care about the denominator for phylogenetic purposes.
    In this function we reformulate the denominator.
    """
    # Mathematica notation:
    # Exp[-c*(1+d)^2 / (2*d)] / (2*d) ) * (
    #       (1+d)*Hypergeometric1F1[1/2, 3/2, c*(1+d)^2 / (2*d)] -
    #       (1-d)*Hypergeometric1F1[1/2, 3/2, c*(1-d)^2 / (2*d)] )
    #
    # precompute a common factor
    c2d = c / (2.*d)
    #
    # compute the asymmetric part
    asym_part = algopy.exp(-c)
    #
    # compute the part that is symmetric in the sense that f(c,d) = f(-c,-d).
    sym_a = 1. / (2.*d)
    sym_b = algopy.exp(-c2d*(d*d + 1.))
    hyper_a = (1. + d) * algopy.special.hyp1f1(0.5, 1.5, c2d*(1+d)**2)
    hyper_b = (1. - d) * algopy.special.hyp1f1(0.5, 1.5, c2d*(1-d)**2)
    sym_part = sym_a * sym_b * (hyper_a - hyper_b)
    #
    # return the approximate value of the function
    return asym_part * sym_part

def kimura_1957_54_denominator_analytic_asymptotic(c, d):
    """
    The numerator of (5.4) goes to 1 in the scaling limit.
    So we only care about the denominator for phylogenetic purposes.
    In this function we reformulate the denominator.
    This uses a large-|z| asymptotic series expansion of hyp1f1.
    """
    a0 = 1. / (2.*c)
    b01 = 1. / (1. + d)
    b02, e = scipy.special.hyp2f0(1.0, 0.5, (2.*d)/(c*(1.+d)*(1.+d)), 2)
    b11 = math.exp(-2.*c) / (1. - d)
    b12, e = scipy.special.hyp2f0(1.0, 0.5, (2.*d)/(c*(1.-d)*(1.-d)), 2)
    return a0 * (b01 * b02 - b11 * b12)

def kimura_genic(c):
    return (1 - math.exp(-2*c)) / (2*c)

def get_form():
    """
    @return: the body of a form
    """
    return [
            ]

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):
    #numpy.set_printoptions(linewidth=200)
    out = StringIO()
    #
    print >> out, 'generating random selection and dominance...'
    print >> out
    for i in range(10):
        #c = random.random()
        #d = random.random()
        c = random.gauss(0, 1)
        d = random.gauss(0, 0.01)
        print >> out, 'c:', c
        print >> out, 'd:', d
        print >> out
        print >> out, 'method a:'
        print >> out, kimura_1957_54_denominator_analytic(c, d)
        print >> out
        print >> out, 'method b:'
        print >> out, kimura_1957_54_denominator_analytic_b(c, d)
        print >> out
        print >> out, 'asymptotic:'
        print >> out, kimura_1957_54_denominator_analytic_asymptotic(c, d)
        print >> out
        print >> out, 'asymptotic, genic:'
        print >> out, kimura_1957_54_denominator_analytic_asymptotic(c, 0)
        print >> out
        print >> out, 'explicitly genic:'
        print >> out, kimura_genic(c)
        print >> out
        print >> out
    #
    return out.getvalue()

