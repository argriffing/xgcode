"""
This module implements functions from Kimura 1957.

It gives explicit analytic solutions of some definite integrals.
Some of the explicit solutions may not be as well behaved numerically as
the numerical integration provided by scipy.integrate.quad.
"""

import unittest

import numpy
import scipy
import scipy.integrate
import mpmath
import algopy
import algopy.special


###########################################################################
# These functions are for the analytical solution of a definite integral.

def denom_not_genic(c, d):
    if not d:
        return numpy.nan
    c2d = c / (2.*d)
    asym_part = algopy.exp(-c)
    sym_a = 1. / (2.*d)
    sym_b = algopy.exp(-c2d*(d*d + 1.))
    hyper_a = (1. + d) * algopy.special.dpm_hyp1f1(0.5, 1.5, c2d*(1+d)**2)
    hyper_b = (1. - d) * algopy.special.dpm_hyp1f1(0.5, 1.5, c2d*(1-d)**2)
    sym_part = sym_a * sym_b * (hyper_a - hyper_b)
    return asym_part * sym_part

def denom_near_genic(c, d):
    """
    This function is better when both |d|<<1 and |d/c|<<1.
    """
    if not c:
        return numpy.nan
    a0 = 1. / (2.*c)
    b01 = 1. / (1.+d)
    b02 = algopy.special.dpm_hyp2f0(1.0, 0.5, (2.*d)/(c*(1.+d)**2))
    b11 = algopy.exp(-2.*c) / (1.-d)
    b12 = algopy.special.dpm_hyp2f0(1.0, 0.5, (2.*d)/(c*(1.-d)**2))
    return a0 * (b01 * b02 - b11 * b12)

def denom_genic_a(c):
    return algopy.special.dpm_hyp1f1(1., 2., -2.*c)

def denom_genic_b(c):
    return (1. - algopy.exp(-2*c)) / (2*c)

def denom_neutral():
    return 1.

def denom_piecewise(c, d):
    """
    This glues together the analytical solution.
    It seems to be usually the case that either denom_near_genic
    or denom_not_genic will give a good solution to the integral,
    but I have not yet found a good criterion for switching between them.
    This is a second attempt, and this time it is
    with respect to the mpmath hypergeometric function implementations.
    """
    eps = 1e-8
    if abs(c) < eps:
        return denom_neutral()
    elif abs(d) < eps:
        return denom_genic_a(c)
    elif -1 < d/c < 1e-3:
        return denom_near_genic(c, d)
    else:
        return denom_not_genic(c, d)


###########################################################################
# These functions are for the numerical solution of a definite integral.


def denom_integrand(x, c, d):
    return algopy.exp(-2*c*d*x*(1-x) - 2*c*x)

def denom_quad(c, d):
    result = scipy.integrate.quad(
            denom_integrand,
            0., 1.,
            args=(c,d),
            full_output=1,
            )
    return result[0]


def test_analytic_integration_solution():
    for c in numpy.linspace(-3, 3, 11):
        for d in numpy.linspace(-0.05, 0.05, 21):
            x = denom_piecewise(c, d)
            y = denom_quad(c, d)
            z = d**2 + (d/c)**2
            print 'c:         ', c
            print 'd:         ', d
            print 'quad:      ', y
            print 'piecewise: ', x
            print 'method:    ', z
            print denom_not_genic(c, d)
            print denom_near_genic(c, d)
            if abs(y - x) / y < 1e-6:
                print 'ok'
            else:
                print '*** bad ***'
            print
    raise Exception

if __name__ == '__main__':
    unittest.main()

