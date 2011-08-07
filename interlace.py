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
import scipy
import sympy
from sympy import matrices
from sympy import abc

import tikz
import sympyutils
import iterutils
import color
import pcurve
import bezier


class Shape:
    """
    This is like a tree or a curve embedded in Euclidean space.
    It is basically a bunch of 1D curves glued together
    and embedded in a higher dimensional Euclidean space.
    Shapes are expected to have the following functions --
    a function that gives bounding box info,
    a function that returns a collection of bezier paths,
    a function that returns the orthoplanar intersections.
    """
    pass

class DifferentiableShape(Shape):
    """
    This is a differentiable parametric curve of interlacing functions.
    """
    def __init__(self, position_exprs, t_initial, t_final):
        """
        The position expressions are univariate sympy expressions.
        Each gives the position along an axis as a function of time t.
        @param position_exprs: sympy expressions that give position at time t
        @param t_initial: initial time
        @param t_final: final time
        """
        velocity_exprs = [x.diff(sympy.abc.t) for x in position_exprs]
        self.fps = [sympyutils.WrappedUniExpr(x) for x in position_exprs]
        self.fvs = [sympyutils.WrappedUniExpr(x) for x in velocity_exprs]
        self.fp = Multiplex(self.fps)
        self.fv = Multiplex(self.fvs)
        self.t_initial = t_initial
        self.t_final = t_final
    def get_bb_min(self):
        """
        Get the min value on each axis.
        """
        values = []
        for f in self.fps:
            t = scipy.optimize.fminbound(f, self.t_initial, self.t_final)
            values.append(f(t))
        return values
    def get_bb_max(self):
        """
        Get the max value on each axis.
        """
        values = []
        for f in self.fps:
            t = scipy.optimize.fminbound(
                    (lambda x: -f(x)), self.t_initial, self.t_final)
            values.append(f(t))
        return values
    def get_orthoplanar_intersections(self):
        """
        Get the list of intersection times on each axis.
        """
        root_seqs = [[]]
        for f in self.fps:
            root_seq = []
            for low, high in iterutils.pairwise(
                    [self.t_initial] + root_seqs[-1] + [self.t_final]):
                root = scipy.optimize.brentq(f, low, high)
                root_seq.append(root)
            root_seqs.append(root_seq)
        return root_seqs[1:]
    def get_bezier_path(self, nchunks=20):
        """
        @param nchunks: use this many chunks in the piecewise approximation
        @return: a BezierPath
        """
        return pcurve.get_bezier_path(
                self.fp, self.fv, self.t_initial, self.t_final, nchunks)


class CubicPolyShape(Shape):
    """
    A parametric cubic polynomial is exactly represented by a Bezier curve.
    Polynomials are sympy Poly objects.
    """
    def __init__(self, polys, t_initial, t_final):
        """
        @param polys: cubic sympy Poly objects
        @param t_initial: initial time
        @param t_final: final time
        """
        self.polys = polys
        self.fps = [sympyutils.WrappedUniPoly(p) for p in polys]
        self.fvs = [sympyutils.WrappedUniPoly(p.diff()) for p in polys]
        self.fp = Multiplex(self.fps)
        self.fv = Multiplex(self.fvs)
        self.t_initial = t_initial
        self.t_final = t_final
    def get_bb_min(self):
        """
        Get the min value on each axis.
        """
        values = []
        for poly in self.polys:
            t, v = sympyutils.poly_fminbound_pair(
                    poly, self.t_initial, self.t_final)
            values.append(v)
        return values
    def get_bb_max(self):
        """
        Get the max value on each axis.
        """
        values = []
        for poly in self.polys:
            t, v = sympyutils.poly_fminbound_pair(
                    -poly, self.t_initial, self.t_final)
            values.append(-v)
        return values
    def get_orthoplanar_intersections(self):
        """
        Get the list of intersection times on each axis.
        """
        return [p.nroots() for p in self.polys]
    def get_bezier_path(self):
        b = bezier.create_bchunk_hermite(
                self.t_initial, self.t_final,
                self.fp(self.t_initial), self.fp(self.t_final),
                self.fv(self.t_initial), self.fv(self.t_final),
                pcurve.OwnedBezierChunk)
        bpath = pcurve.BezierPath([b])
        b.parent_ref = bpath
        return bpath


class DifferentialCubicPolyShape(CubicPolyShape):
    """
    The polynomials defining the shape are successive derivatives.
    """
    def __init__(self, cubic_poly, t_initial, t_final):
        """
        @param cubic_poly: a cubic sympy Poly with distinct roots
        @param t_initial: initial time
        @param t_final: final time
        """
        self.fps = fps
        self.t_initial = t_initial
        self.t_final = t_final

class PiecewiseLinearPathShape(Shape):
    """
    A path of line segments in higher dimensional Euclidean space.
    """
    def get_bb_min(self):
        pass
    def get_bb_max(self):
        pass

class PiecewiseLinearTreeShape(Shape):
    """
    A path of line segments in higher dimensional Euclidean space.
    """
    def get_bb_min(self):
        pass
    def get_bb_max(self):
        pass





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
    """
    Turn a sequence of functions into a single function.
    Each component function should return a float given a float,
    while the returned function will return a numpy array given a float.
    """
    def __init__(self, fs):
        """
        @param fs: sequence of (float -> float) python functions
        """
        self.fs = fs
    def __call__(self, t):
        return np.array([f(t) for f in self.fs])

class MultiplexPolys:
    def __init__(self, polys):
        """
        @param fns: sequence of univariate sympy Poly objects
        """
        self.polys = polys
    def __call__(self, t):
        return np.array([float(p.eval(t)) for p in self.polys])



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

