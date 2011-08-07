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


class TestSympyUtils:

    def test_wrapped_univariate_expression(self):
        expr = sympy.abc.t ** 2.0 + 1.0
        f = WrappedUniExpr(expr)
        observed = f(3.0)
        expected = 10.0
        self.assertEqual(observed, expected)


if __name__ == '__main__':
    unittest.main()

