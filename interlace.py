"""
Interlacing functions.

This includes polynomials.
Some families of polynomial sequences are known to
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

#TODO use sympy for characteristic polynomials and differentiation

class InterlacingPoly:
    def __init__(self, root_a, root_b, root_c, initial_t, final_t):
        """
        @param root_a: first root of monic cubic polynomial
        @param root_b: second root of monic cubic polynomial
        @param root_c: third root of monic cubic polynomial
        """
        self.root_a = root_a
        self.root_b = root_b
        self.root_c = root_c
    def get_cubic_roots(self):
        a = self.root_a
        b = self.root_b
        c = self.root_c
        return a, b, c
    def get_quadratic_roots(self):
        a = self.root_a
        b = self.root_b
        c = self.root_c
        A = a*a + b*b + c*c
        B = a*b + a*c + b*c
        S = a + b + c
        r0 = float(S + math.sqrt(A - B)) / 3
        r1 = float(S - math.sqrt(A - B)) / 3
        return r0, r1
    def get_linear_root(self):
        r = float(self.root_a + self.root_b + self.root_c) / 3
        return r
    def __call__(self, t):
        """
        @param t: a float
        @return: a 3d point
        """
        a = self.root_a
        b = self.root_b
        c = self.root_c
        z = (t - a) * (t - b) * (t - c)
        y = 3*t*t - 2*t*(a + b + c) + a*b + a*c + b*c
        x = 6*t - 2*(a + b + c)
        return np.array([x, y, z])


class TestInterlacing(unittest.TestCase):
    
    def test_default(self):
        M = sympy.matrices.Matrix((
            (1, 2),
            (2, 1)))
        print M.berkowitz()
        print M.berkowitz_charpoly(sympy.abc.x)
        print M.charpoly(sympy.abc.x)


if __name__ == '__main__':
    unittest.main()

