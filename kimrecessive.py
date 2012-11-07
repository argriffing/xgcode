"""
This module implements functions from Kimura 1957.

It gives explicit analytic solutions of some definite integrals.
Some of the explicit solutions may not be as well behaved numerically as
the numerical integration provided by scipy.integrate.quad.
"""

import numpy
from numpy import testing
import scipy
import scipy.integrate
import algopy
import algopy.special



def _hyp2f0_combo(a1, a2, x):
    """
    This is a helper function.
    It works around a bugginess of scipy hyp2f0 implementation.
    """
    if x < -1e-3:
        return (-1/x)**a1 * algopy.special.hyperu(a1, 1+a1-a2, (-1/x))
    else:
        return algopy.special.hyp2f0(a1, a2, x)


###########################################################################
# These functions are for the analytical solution of a definite integral.

def denom_knudsen(c):
    """
    This is +gwF = 1/2.
    """
    return algopy.exp(-c)

def denom_complete_dominant(c):
    return algopy.special.hyp1f1(1.0, 1.5, -2*c)

def denom_complete_recessive(c):
    return algopy.special.hyp1f1(0.5, 1.5, -2*c)

def denom_not_genic(c, d):
    #if not d:
        #return numpy.nan
    c2d = c / (2.*d)
    asym_part = algopy.exp(-c)
    sym_a = 1. / (2.*d)
    sym_b = algopy.exp(-c2d*(d*d + 1.))
    hyper_a = (1. + d) * algopy.special.dpm_hyp1f1(0.5, 1.5, c2d*(1+d)**2)
    hyper_b = (1. - d) * algopy.special.dpm_hyp1f1(0.5, 1.5, c2d*(1-d)**2)
    sym_part = sym_a * sym_b * (hyper_a - hyper_b)
    return asym_part * sym_part

def denom_near_genic_combo(c, d):
    a0 = 1. / (2.*c)
    b01 = 1. / (1.+d)
    z1 = (2.*d)/(c*(1.+d)**2)
    z2 = (2.*d)/(c*(1.-d)**2)
    z1_recip = (c*(1.+d)**2) / (2.*d)
    #FIXME: unfinished
    b02 = _hyp2f0_combo(1.0, 0.5, (2.*d)/(c*(1.+d)**2))
    b11 = algopy.exp(-2.*c) / (1.-d)
    b12 = _hyp2f0_combo(1.0, 0.5, (2.*d)/(c*(1.-d)**2))
    return a0 * (b01 * b02 - b11 * b12)

def denom_near_genic(c, d):
    #if not c:
        #return numpy.nan
    #if d in (-1, 1):
        #return numpy.nan
    a0 = 1. / (2.*c)
    b01 = 1. / (1.+d)
    b02 = algopy.special.dpm_hyp2f0(1.0, 0.5, (2.*d)/(c*(1.+d)**2))
    b11 = algopy.exp(-2.*c) / (1.-d)
    b12 = algopy.special.dpm_hyp2f0(1.0, 0.5, (2.*d)/(c*(1.-d)**2))
    return a0 * (b01 * b02 - b11 * b12)

def denom_genic_a(c):
    return algopy.special.hyp1f1(1., 2., -2.*c)

def denom_genic_b(c):
    return (1. - algopy.exp(-2*c)) / (2*c)

def denom_neutral():
    return 1.

def denom_piecewise(c, d):
    """
    This glues together the analytical solution.
    This is a second attempt, and this time it is
    with respect to the mpmath hypergeometric function implementations.
    """
    small_eps = 1e-8
    large_eps = 1e-3
    if abs(c) < small_eps:
        return denom_neutral()
    elif abs(d) < small_eps:
        return denom_genic_a(c)
    elif abs(d) > 1 - large_eps:
        return denom_not_genic(c, d)
    elif -1 < d/c < large_eps:
        return denom_near_genic(c, d)
    else:
        return denom_not_genic(c, d)

def denom_hyperu_b(c, d):
    """
    This uses only algopy.special.hyperu(1.0, 1.5, x).
    It is only meaningful when c and d have different signs,
    but for now do not enforce this.
    """
    prefix = algopy.exp(-c)
    a1 = algopy.exp(-c)
    a2 = (1. - d) / (4. * d)
    a3 = algopy.special.hyperu(1.0, 1.5, -c*(1-d)**2 / (2. * d))
    b1 = algopy.exp(c)
    b2 = (1. + d) / (4. * d)
    b3 = algopy.special.hyperu(1.0, 1.5, -c*(1+d)**2 / (2. * d))
    return prefix * (a1 * a2 * a3 - b1 * b2 * b3)

def denom_hyp1f1_b(c, d):
    """
    This uses only algopy.special.hyp1f1(1.0, 1.5, x).
    """
    prefix = algopy.exp(-c)
    a1 = algopy.exp(c)
    a2 = (1. + d) / (2. * d)
    a3 = algopy.special.hyp1f1(1.0, 1.5, -c*(1+d)**2 / (2. * d))
    b1 = algopy.exp(-c)
    b2 = (1. - d) / (2. * d)
    b3 = algopy.special.hyp1f1(1.0, 1.5, -c*(1-d)**2 / (2. * d))
    return prefix * (a1 * a2 * a3 - b1 * b2 * b3)

def denom_hyp2f0_b(c, d):
    prefix = algopy.exp(-c)
    a1 = algopy.exp(c)
    a2 = 1. / (2. * c)
    a3 = 1. / (1. + d)
    a4 = algopy.special.hyp2f0(1.0, 0.5, (2. * d) / (c * (1. + d)**2))
    b1 = algopy.exp(-c)
    b2 = 1. / (2. * c)
    b3 = 1. / (1. - d)
    b4 = algopy.special.hyp2f0(1.0, 0.5, (2. * d) / (c * (1. - d)**2))
    return prefix * (a1 * a2 * a3 * a4 - b1 * b2 * b3 * b4)

def denom_dawsn_b(c, d):
    p1 = scipy.exp(-c)
    p2 = 1. / (scipy.sqrt(c) * scipy.sqrt(2*d))
    a1 = scipy.exp(-c)
    a2 = scipy.special.dawsn((d-1)*scipy.sqrt(c) / scipy.sqrt(2*d))
    b1 = scipy.exp(c)
    b2 = scipy.special.dawsn((d+1)*scipy.sqrt(c) / scipy.sqrt(2*d))
    return p1 * p2 * (a1 * a2 + b1 * b2)

def denom_erfi_b(c, d):
    p1 = scipy.sqrt(scipy.pi)
    p2 = 1. / (2 * scipy.sqrt(c) * scipy.sqrt(2*d))
    p3 = scipy.exp(-(c*(1.+d)*(1.+d)) / (2.*d))
    a = scipy.special.erfi((scipy.sqrt(c)*(1.+d))/scipy.sqrt(2*d))
    b = scipy.special.erfi((scipy.sqrt(c)*(1.-d))/scipy.sqrt(2*d))
    return p1 * p2 * p3 * (a - b)


###########################################################################
# These functions are for the numerical solution of a definite integral.


def denom_integrand(x, c, d):
    return numpy.exp(-2*c*d*x*(1-x) - 2*c*x)

def denom_quad(c, d):
    result = scipy.integrate.quad(
            denom_integrand,
            0., 1.,
            args=(c,d),
            full_output=1,
            )
    return result[0]

class Test_KimuraRecessive(testing.TestCase):

    def test_neutral(self):
        c = 0.0
        for d in (-0.123, 0.01, 1.23):
            x = denom_neutral()
            y = denom_genic_a(c)
            z = denom_not_genic(c, d)
            w = denom_quad(c, d)
            testing.assert_allclose(x, y)
            testing.assert_allclose(x, z)
            testing.assert_allclose(x, w)

    def test_genic(self):
        d = 0.0
        for c in (-0.123, 0.01, 1.23):
            x = denom_genic_a(c)
            y = denom_genic_b(c)
            z = denom_near_genic(c, d)
            w = denom_quad(c, d)
            testing.assert_allclose(x, y)
            testing.assert_allclose(x, z)
            testing.assert_allclose(x, w)

    def test_complete_dominant(self):
        d = 1.0
        for c in (-0.123, 0.01, 1.23):
            x = denom_complete_dominant(c)
            y = denom_not_genic(c, d)
            w = denom_quad(c, d)
            testing.assert_allclose(x, y)
            testing.assert_allclose(x, w)

    def test_complete_recessive(self):
        d = -1.0
        for c in (-0.123, 0.01, 1.23):
            x = denom_complete_recessive(c)
            y = denom_not_genic(c, d)
            w = denom_quad(c, d)
            testing.assert_allclose(x, y)
            testing.assert_allclose(x, w)

if __name__ == '__main__':
    testing.run_module_suite()

