"""
This is an interface between sympy expressions and python functions.
"""

import unittest

import sympy
from sympy import abc

class WrappedUniExpr:
    """
    This is a wrapped univariate sympy expression.
    The idea is to make this expression available as a python function.
    The variable is t.
    """
    def __init__(self, sympy_expr):
        self.sympy_expr = sympy_expr
    def __call__(self, t):
        """
        @param t: a python float
        @return: a python float
        """
        return float(self.sympy_expr.subs(sympy.abc.t, t))

class WrappedUniPoly:
    """
    This is a wrapped univariate sympy Poly object.
    The idea is to make this polynomial available as a python function.
    The variable is t.
    """
    def __init__(self, sympy_poly_object):
        self.sympy_poly_object = sympy_poly_object
    def __call__(self, t):
        """
        @param t: a python float
        @return: a python float
        """
        return float(self.sympy_poly_object.eval(t))

def roots_to_poly(roots):
    """
    @param roots: a collection of real roots
    @return: a sympy Poly object with parameter t
    """
    t = sympy.abc.t
    p = sympy.Poly(1.0, t)
    for r in roots:
        p *= (t - r)
    return p

def poly_fminbound_pair(p, low, high):
    """
    Return the (time, value) pair.
    For max bound use -fminbound(-p).
    @param p: a sympy Poly object
    @param low: low end of the interval
    @param high: high end of the interval
    @return: the (t, v) pair where v is minimum
    """
    vt_pairs = []
    for t in (low, high):
        vt_pairs.append((float(p.eval(t)), t))
    for t in p.diff().nroots():
        vt_pairs.append((float(p.eval(t)), t))
    v, t = min(vt_pairs)
    return (t, v)

def poly_fminbound(p, low, high):
    """
    Return the time only.
    For max bound use -fminbound(-p).
    @param p: a sympy Poly object
    @param low: low end of the interval
    @param high: high end of the interval
    @return: the time t at which the min is attained
    """
    t, v = poly_fminbound_pair(p, low, high)
    return t


class TestSympyUtils(unittest.TestCase):

    def test_wrapped_univariate_expression(self):
        expr = sympy.abc.t ** 2 + 1.0
        f = WrappedUniExpr(expr)
        observed = f(3.0)
        expected = 10.0
        self.assertEqual(observed, expected)

    def test_wrapped_univariate_polynomial_object(self):
        poly = sympy.Poly(sympy.abc.t ** 2 + 1.0, sympy.abc.t)
        f = WrappedUniPoly(poly)
        observed = f(3.0)
        expected = 10.0
        self.assertEqual(observed, expected)

    def make_sample_cubic(self):
        """
        The cubic poly has local max and 1 and local min at 2.
        """
        ten_thirds = sympy.Rational(10,3)
        t = sympy.abc.t
        return sympy.Poly(
                ten_thirds * t**3 - 15*t**2 + 20*t - 7, t, domain='RR')

    def test_poly_fminbound_low(self):
        p = self.make_sample_cubic()
        low = 0.1
        high = 2.1
        expected = 0.1
        observed = poly_fminbound(p, low, high)
        self.assertAlmostEqual(observed, expected)

    def test_poly_fminbound_high(self):
        p = -self.make_sample_cubic()
        low = 0.9
        high = 3.0
        expected = 3.0
        observed = poly_fminbound(p, low, high)
        self.assertAlmostEqual(observed, expected)

    def test_poly_fminbound_local_min(self):
        p = self.make_sample_cubic()
        low = 0.9
        high = 2.1
        expected = 2.0
        observed = poly_fminbound(p, low, high)
        self.assertAlmostEqual(observed, expected)

    def test_principal_submatrix_charpoly(self):
        t = sympy.abc.t
        M = sympy.Matrix([
            [5, -4, 0],
            [-4, 12, -5],
            [0, -5, 7]])
        tI = t*sympy.eye(3)
        p3x3 = (tI - M).det()
        p2x2 = (tI - M)[:2, :2].det()
        p1x1 = (tI - M)[:1, :1].det()
        p3x3_expected = t**3 - 24*(t**2) + 138*t - 183
        p2x2_expected = t**2 - 17*t + 44
        p1x1_expected = t - 5
        self.assertEqual(p3x3, p3x3_expected)
        self.assertEqual(p2x2, p2x2_expected)
        self.assertEqual(p1x1, p1x1_expected)

    def test_schur_complement_charpoly(self):
        t = sympy.abc.t
        M = sympy.Matrix([
            [5, -4, 0],
            [-4, 12, -5],
            [0, -5, 7]])
        tI = t*sympy.eye(3)
        comp2x2 = (M[:2, :2] - M[:2, 2]*(1/M[2,2])*M[2, :2])
        comp1x1 = (M[:1, :1] - M[0, 1:]*M[1:, 1:].inv()*M[1:, 0])
        p3x3 = (tI - M).det()
        p2x2 = (tI[:2, :2] - comp2x2).det()
        p1x1 = (tI[:1, :1] - comp1x1).det()
        p3x3_expected = t**3 - 24*(t**2) + 138*t - 183
        p2x2_expected = t**2 - t*sympy.Rational(94, 7) + sympy.Rational(183, 7)
        p1x1_expected = t - sympy.Rational(183, 59)
        self.assertEqual(p3x3, p3x3_expected)
        self.assertEqual(p2x2, p2x2_expected)
        self.assertEqual(p1x1, p1x1_expected)


if __name__ == '__main__':
    unittest.main()

