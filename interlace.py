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
import tikz
import iterutils
import color

def is_strictly_increasing(seq):
    for a, b in iterutils.pairwise(seq):
        if not a < b:
            return False
    return True

def assert_support(t_seq, y_seqs):
    # check that the t sequence is increasing
    if not is_strictly_increasing(t_seq):
        msg = 'expected a strictly increasing t sequence'
        raise ValueError(msg)
    # check that each sequence has the same number of samples
    seqs = [t_seq] + y_seqs
    lengths = set(len(seq) for seq in seqs)
    if len(lengths) != 1:
        msg = 'expected each sequence to have the same length'
        raise ValueError(msg)

def tikz_superposition(t_seq, y_seqs, width, height):
    """
    Return the body of a tikzpicture environment.
    The input defines k piecewise-linear parametric functions.
    The domain is a real interval.
    The kth sequence of y values should have k zero-crossings.
    The returned drawing has arbitrary horizontal and vertical scale.
    @param t_seq: sequence of t values
    @param y_seqs: sequence of y value sequences
    @param width: a horizontal scaling factor
    @param height: a vertical scaling factor
    @return: LaTeX code for a tikzpicture
    """
    # check the form of the input sequences
    assert_support(t_seq, y_seqs)
    # Get the y scaling factor.
    # This is the value by which the observed y values are multiplied
    # to give a number that is relevant to the tikz coordinate system.
    ymin = min(min(seq) for seq in y_seqs)
    ymax = max(max(seq) for seq in y_seqs)
    yscale = height / float(ymax - ymin)
    # Get the x scaling factor.
    xmin = t_seq[0]
    xmax = t_seq[-1]
    xscale = width / float(xmax - xmin)
    # Get the rescaled sequences.
    t_seq_rescaled = [t*xscale for t in t_seq]
    y_seqs_rescaled = [[y*yscale for y in seq] for seq in y_seqs]
    # Start the text array.
    arr = []
    # Plot the horizontal domain segment.
    pa = tikz.point_to_tikz((t_seq_rescaled[0], 0))
    pb = tikz.point_to_tikz((t_seq_rescaled[-1], 0))
    arr.append('\\draw %s -- %s;' % (pa, pb))
    # Plot the scaled functions.
    for seq, c in zip(y_seqs_rescaled, color.wolfram):
        arr.append('\\draw[color=%s]' % c)
        points = zip(t_seq_rescaled, seq)
        arr.append(tikz.curve_to_tikz(points, 3) + ';')
    # Return the text.
    return '\n'.join(arr)


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
        #return np.array([float(f(t)) for f in self.fns])
        return np.array([float(f.eval(t)) for f in self.fns])


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

    def test_ndiff_roots_symbolic(self):
        roots = (1.25, 4.5, 5.75)
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

