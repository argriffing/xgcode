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
from scipy import linalg
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
    def get_infinity_radius(self):
        """
        Infinity refers to the infinity norm.
        @return: max of absolute values of axis aligned bounding box coords
        """
        a = np.linalg.norm(self.get_bb_min(), np.Inf)
        b = np.linalg.norm(self.get_bb_max(), np.Inf)
        return max(a, b)

class ParametricShape(Shape):
    """
    This 1D shape is parameterized by a single variable.
    The fp member function should map the time to a point.
    """
    def get_bb_min(self):
        """
        Get the min value on each axis.
        """
        times = self.get_bb_min_times()
        points = [self.fp(t) for t in times]
        return np.min(points, axis=0)
    def get_bb_max(self):
        """
        Get the max value on each axis.
        """
        times = self.get_bb_max_times()
        points = [self.fp(t) for t in times]
        return np.max(points, axis=0)
    def get_orthoplanar_intersections(self):
        """
        Get the list of intersection points per axis.
        """
        point_seqs = []
        for time_seq in self.get_orthoplanar_intersection_times():
            point_seqs.append([self.fp(t) for t in time_seq])
        return point_seqs
    def get_bezier_paths(self):
        return [self.get_bezier_path()]

class DifferentiableShape(ParametricShape):
    """
    This is a differentiable parametric curve of interlacing functions.
    """
    def __init__(self, position_exprs, t_initial, t_final, nchunks_default=20):
        """
        The position expressions are univariate sympy expressions.
        Each gives the position along an axis as a function of time t.
        @param position_exprs: sympy expressions that give position at time t
        @param t_initial: initial time
        @param t_final: final time
        @param nchunks_default: the default number of chunks for a bpath
        """
        velocity_exprs = [x.diff(sympy.abc.t) for x in position_exprs]
        self.fps = [sympyutils.WrappedUniExpr(x) for x in position_exprs]
        self.fvs = [sympyutils.WrappedUniExpr(x) for x in velocity_exprs]
        self.fp = Multiplex(self.fps)
        self.fv = Multiplex(self.fvs)
        self.t_initial = t_initial
        self.t_final = t_final
        self.nchunks_default = nchunks_default
    def get_bb_min_times(self):
        """
        Get the time of the min value on each axis.
        """
        return [scipy.optimize.fminbound(
            f, self.t_initial, self.t_final) for f in self.fps]
    def get_bb_max_times(self):
        """
        Get the time of the max value on each axis.
        """
        return [scipy.optimize.fminbound(
            (lambda x: -f(x)), self.t_initial, self.t_final) for f in self.fps]
    def get_orthoplanar_intersection_times(self):
        """
        Get the intersection times for the plane orthogonal to each axis.
        Note that this function assumes interlacing roots.
        """
        root_seqs = [[]]
        for f in self.fps:
            root_seq = []
            for low, high in iterutils.pairwise(
                    [self.t_initial] + root_seqs[-1] + [self.t_final]):
                root_seq.append(scipy.optimize.brentq(f, low, high))
            root_seqs.append(root_seq)
        return root_seqs[1:]
    def get_bezier_path(self, nchunks_in=None):
        """
        @param nchunks_in: use this many chunks in the piecewise approximation
        @return: a BezierPath
        """
        nchunks = self.nchunks_default
        if nchunks_in:
            nchunks = nchunks_in
        return pcurve.get_bezier_path(
                self.fp, self.fv, self.t_initial, self.t_final, nchunks)

class CubicPolyShape(ParametricShape):
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
    def get_bb_min_times(self):
        """
        Get the time of the min value on each axis.
        """
        return [sympyutils.poly_fminbound(
            poly, self.t_initial, self.t_final) for poly in self.polys]
    def get_bb_max_times(self):
        """
        Get the time of the max value on each axis.
        """
        return [sympyutils.poly_fminbound(
            -poly, self.t_initial, self.t_final) for poly in self.polys]
    def get_orthoplanar_intersection_times(self):
        """
        Get the intersection times for the plane orthogonal to each axis.
        """
        return [p.nroots() for p in self.polys]
    def get_bezier_path(self):
        b = bezier.create_bchunk_hermite(
                self.t_initial, self.t_final,
                self.fp(self.t_initial), self.fp(self.t_final),
                self.fv(self.t_initial), self.fv(self.t_final))
        return pcurve.BezierPath([b])

class PiecewiseLinearPathShape(Shape):
    def __init__(self, points):
        """
        @param points: a sequence of high dimensional points as numpy arrays
        """
        self.points = [np.array(p) for p in points]
        self.ndim = len(points[0])
    def get_bb_min(self):
        """
        Get the min value on each axis.
        """
        return np.min(self.points, axis=0)
    def get_bb_max(self):
        """
        Get the max value on each axis.
        """
        return np.max(self.points, axis=0)
    def get_orthoplanar_intersections(self):
        """
        Get the list of intersection points per axis.
        """
        abstol = 1e-6
        point_seqs = []
        for axis in range(self.ndim):
            point_seq = []
            # check points for exact intersections
            for p in self.points:
                if abs(p[axis]) < abstol:
                    point_seq.append(p)
            # check line segments for intersections
            for pa, pb in iterutils.pairwise(self.points):
                if abs(pa[axis]) > abstol and abs(pb[axis]) > abstol:
                    if pa[axis]*pb[axis] < 0:
                        p = (pb[axis]*pa - pa[axis]*pb) / (pb[axis] - pa[axis])
                        point_seq.append(p)
            point_seqs.append(point_seq)
        return point_seqs
    def get_bezier_path(self):
        bchunks = []
        for i, (pa, pb) in enumerate(iterutils.pairwise(self.points)):
            b = bezier.create_bchunk_line_segment(pa, pb)
            b.start_time = float(i)
            b.stop_time = float(i+1)
            bchunks.append(b)
        return pcurve.BezierPath(bchunks)
    def get_bezier_paths(self):
        return [self.get_bezier_path()]

class PiecewiseLinearTreeShape(Shape):
    def __init__(self, wat):
        pass


class Hypershape:
    """
    A hypershape defines points in a more-than-geometric space.
    For example each point could be like (x1, x2, y1, y2, y3) meaning
    that at (x1, x2) the value is (y1, y2, y3).
    A purely geometric shape would only care about
    the set of (y1, y2, y3) triples and would not care about the domain.
    """
    pass

class PiecewiseLinearPathHypershape(Hypershape):
    """
    This is a sequence of points augmented with timing information.
    Practically, this should be used for the superposition figures and the 
    sign cut figures.
    """
    pass

class PiecewiseLinearTreeHypershape(Hypershape):
    """
    This is a set of points with 2d tree topology and layout information.
    Practically, this should be used for the superposition figures and the 
    sign cut figures.
    """
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

def get_tikz_bezier_2d(bpath):
    lines = []
    # draw everything except for the last point of the last chunk
    for b in bpath.bchunks:
        pts = [tikz.point_to_tikz(p) for p in b.get_points()[:-1]]
        lines.append('%s .. controls %s and %s ..' % tuple(pts))
    # draw the last point of the last chunk
    lines.append('%s;' % tikz.point_to_tikz(bpath.bchunks[-1].p3))
    return '\n'.join(lines)

def tikz_shape_superposition(shapes, width, height):
    """
    Return the body of a tikzpicture environment.
    @param shapes: a sequence of 2D Shape objects
    @param width: max tikz width
    @param height: max tikz height
    @return: tikz text
    """
    bbmax = np.max([shape.get_bb_max() for shape in shapes], axis=0)
    bbmin = np.min([shape.get_bb_min() for shape in shapes], axis=0)
    scale = np.array([width, height], dtype=float) / (bbmax - bbmin)
    f = lambda x: x*scale
    colors = ['black'] + color.wolfram
    arr = []
    for c, shape in zip(colors, shapes):
        for bpath in shape.get_bezier_paths():
            bpath.transform(f)
            # Chop up the bpath so that the stupid bounding box that tikz
            # makes will still be approximately correct.
            bpath.refine_for_bb()
            # Draw the bpath.
            arr.extend([
                '\\draw[thick,%s]' % c,
                get_tikz_bezier_2d(bpath)])
    return '\n'.join(arr)

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


def roots_to_differential_polys(roots):
    """
    Construct a sequence of interlacing polynomials.
    The input is the collection of distinct roots
    of the highest degree polynomial in the sequence.
    @param roots: a collection of distinct roots
    @return: a sequence of interlacing polynomials
    """
    sympy_t = sympy.abc.t
    if len(roots) != len(set(roots)):
        raise ValueError('expected distinct roots')
    p = sympyutils.roots_to_poly(roots)
    nroots = len(roots)
    polys = [p]
    for i in range(nroots-1):
        p = polys[-1].diff(sympy_t)
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
    Turn an iterable of functions into a single function.
    Each component function should return a float given a float,
    while the returned function will return a numpy array given a float.
    """
    def __init__(self, fs):
        """
        @param fs: iterable of (float -> float) python functions
        """
        self.fs = list(fs)
    def __call__(self, t):
        return np.array([f(t) for f in self.fs])





class TestInterlacing(unittest.TestCase):

    def test_roots_to_poly(self):
        roots = (1.0, 4.0, 5.0)
        p = sympyutils.roots_to_poly(roots)
        self.assertTrue(p.is_monic)
        root_to_count = sympy.roots(p)
        self.assertEqual(set(root_to_count.values()), set([1]))
        observed = set(root_to_count.keys())
        expected = set(roots)
        for r_observed, r_expected in zip(sorted(observed), sorted(expected)):
            self.assertAlmostEqual(r_observed, r_expected)

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
        polys = roots_to_differential_polys(roots)
        f = Multiplex((sympyutils.WrappedUniPoly(p) for p in polys))
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

