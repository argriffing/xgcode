"""
Interlacing functions.

This includes polynomials.  Some families of polynomial sequences are known to
strictly or weakly interlace.
This includes all sequences of polynomials that you get by
starting with a high degree polynomial with distinct real roots
and then taking derivatives.
Also all orthogonal polynomial sequences interlace.
So do certain sequences of characteristic polynomials.
This includes characteristic polynomials of sequences
of real symmetric matrices whose dimension is reduced by
taking principal submatrices.
Also reducing the dimensions of real symmetric positive semidefinite matrices
by taking the Schur complement gives a sequence
of interlacing characteristic polynomials.
"""

import math
import unittest

import numpy as np
import sympy
from sympy import matrices
from sympy import abc

def roots_to_poly(roots):
    """
    Given some roots, return a monic polynomial.
    @return: a sympy polynomial
    """
    p_prod = sympy.Poly(1, sympy.abc.x)
    for r in roots:
        term = sympy.Poly(sympy.abc.x - r, sympy.abc.x)
        p_prod = p_prod * term
    return p_prod

def roots_to_differential_polys(roots):
    """
    Construct a sequence of interlacing polynomials.
    The input is the collection of distinct roots
    of the highest degree polynomial in the sequence.
    @param roots: a collection of distinct roots
    @return: a sequence of interlacing polynomials
    """
    if len(roots) != len(set(roots)):
        raise ValueError('expected distinct roots')
    p = roots_to_poly(roots)
    nroots = len(roots)
    polys = [p]
    for i in range(nroots-1):
        p = polys[-1].diff(sympy.abc.x)
        polys.append(p)
    polys.reverse()
    return polys

def matrix_to_principal_polys(M):
    """
    Construct a sequence of interlacing polynomials.
    @param M: real symmetric matrix
    @return: a sequence of interlacing polynomials
    """
    pass

def matrix_to_schur_polys(M):
    """
    Construct a sequence of interlacing polynomials.
    @param M: positive semidefinite real symmetric matrix
    @return: a sequence of interlacing polynomials
    """
    pass


class Multiplex:
    def __init__(self, fns):
        """
        @param fns: sequence of univariate functions
        """
        self.fns = fns
    def __call__(self, t):
        return np.array([float(f(t)) for f in self.fns])


class TestInterlacing(unittest.TestCase):

    def test_roots_to_poly(self):
        roots = (1, 4, 5)
        p = roots_to_poly(roots)
        self.assertTrue(p.is_monic)
        root_to_count = sympy.roots(p)
        self.assertEqual(set(root_to_count.values()), set([1]))
        expected = set(sympy.Integer(x) for x in roots)
        self.assertEqual(set(root_to_count.keys()), expected)

    def test_ndiff_roots(self):
        roots = (1, 4, 5)
        a, b, c = roots
        polys = roots_to_differential_polys(roots)
        # compute linear root manually
        r = float(a + b + c) / 3
        # check linear root
        observed = sorted(float(r) for r in sympy.roots(polys[0]))
        expected = sorted([r])
        self.assertTrue(np.allclose(observed, expected))
        # compute quadratic roots manually
        A = a*a + b*b + c*c
        B = a*b + a*c + b*c
        S = a + b + c
        r0 = float(S + math.sqrt(A - B)) / 3
        r1 = float(S - math.sqrt(A - B)) / 3
        # check quadratic roots
        observed = sorted(float(r) for r in sympy.roots(polys[1]))
        expected = sorted([r0, r1])
        self.assertTrue(np.allclose(observed, expected))
        # check cubic roots
        observed = sorted(float(r) for r in sympy.roots(polys[2]))
        expected = sorted(roots)
        self.assertTrue(np.allclose(observed, expected))
    
    def test_poly_eval(self):
        roots = (1, 4, 5)
        t = 3.2
        a, b, c = roots
        f = Multiplex(roots_to_differential_polys(roots))
        # compute the evaluation manually
        z = (t - a) * (t - b) * (t - c)
        y = 3*t*t - 2*t*(a + b + c) + a*b + a*c + b*c
        x = 6*t - 2*(a + b + c)
        # check the evaluation
        observed = f(t)
        expected = np.array([x, y, z])
        self.assertTrue(np.allclose(observed, expected))



if __name__ == '__main__':
    unittest.main()

