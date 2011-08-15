"""
Vignettes for interlacing.

Most of the hardcoded values in this module should eventually
be replaced by automatically obtained values.
-
Each sample should provide enough information to make four visualizations.
1) a 3d thumbnail logo with a fixed rotation and translation
2) a 3d image with a fixed rotation but variable translation
3) a superposition image (2d for linear shapes, 3d for tree shapes)
4) a four-vertical-panel sign interlacing plot
-
Hardcoded items shared between the thumbnail logo and the 3d image.
1) the interlace.Shape object
-
Hardcoded items for the thumbnail logo.
1) a scaling factor
-
Hardcoded items for the 3d image.
1) a scaling factor
2) six axis radii
-
"""

import math

import numpy as np
import scipy
from scipy import linalg
import sympy
from sympy import matrices
from sympy import abc

import sympyutils
import interlace
import Newick
import FtreeIO
import Ftree

def get_samples():
    return [
            DerivativePoly(),
            PrincipalCharpoly(),
            SchurCharpoly(),
            OrthogonalPoly(),
            SturmLiouville(),
            FiniteDifferences(),
            LaplacePath(),
            LaplaceTree(),
            SchurTree()]

def _get_cubic_superposition_shapes(
        polys, initial_t, final_t):
    """
    This is a helper function.
    """
    shapes = []
    # add the axis shape
    x_axis = interlace.PiecewiseLinearPathShape([
        (initial_t, 0),
        (final_t, 0)])
    shapes.append(x_axis)
    # each 2d shape is a parametric cubic polynomial
    t = sympy.abc.t
    axis_poly = sympy.Poly(t, t)
    for poly in polys:
        shape = interlace.CubicPolyShape(
                (axis_poly, poly), initial_t, final_t)
        shapes.append(shape)
    return shapes

def _get_differentiable_superposition_shapes(
        exprs, initial_t, final_t, nsegs=10):
    """
    This is a helper function.
    @param exprs: sympy expressions
    """
    shapes = []
    # add the axis shape
    x_axis = interlace.PiecewiseLinearPathShape([
        (initial_t, 0),
        (final_t, 0)])
    shapes.append(x_axis)
    # each 2d shape is a differentiable shape
    for expr in exprs:
        shape = interlace.DifferentiableShape(
                (sympy.abc.t, expr), initial_t, final_t, nsegs)
        shapes.append(shape)
    return shapes

class Sample:
    """
    Default values to use for debugging.
    """
    def get_small_3d_sf(self):
        """
        @return: scaling factor for the small 3d image
        """
        return 1.0
    def get_large_3d_sf(self):
        """
        @return: scaling factor for the large 3d image
        """
        return 1.0
    def get_axis_radii(self):
        """
        Return the six half axis radii in tikz units.
        @return: x+ rad, x- rad, y+ rad, y- rad, z+ rad, z- rad
        """
        return [1.0]*6


class DerivativePoly(Sample):
    """
    A parametric curve defined by a cubic polynomial and its derivatives.
    """
    def __init__(self):
        # this is basically like 
        # f(x) = (x-2)(x-5)(x-7)
        # f(x) = x^3 - 14 x^2 + 59 x - 70
        # f'(x) = 3 x^2 - 28 x + 59
        # f''(x) = 6 x - 28
        self.initial_t = 0.9
        root_a = 1.0
        root_b = 2.5
        root_c = 3.5
        self.final_t = 3.6
        p3 = sympyutils.roots_to_poly((root_a, root_b, root_c))
        p2 = p3.diff()
        p1 = p2.diff()
        self.polys = (p1, p2, p3)
        self.shape =  interlace.CubicPolyShape(
                self.polys, self.initial_t, self.final_t)
    def get_shape(self):
        return self.shape
    def get_small_3d_sf(self):
        r = self.shape.get_infinity_radius()
        return 0.8 * (1.0 / r)
    def get_large_3d_sf(self):
        return 1.0
    def get_axis_radii(self):
        xp_rad = 6.0
        xn_rad = 6.0
        yp_rad = 6.0
        yn_rad = 3.0
        zp_rad = 3.0
        zn_rad = 3.0
        return xp_rad, xn_rad, yp_rad, yn_rad, zp_rad, zn_rad
    def get_superposition_shapes(self):
        return _get_cubic_superposition_shapes(
                self.polys, self.initial_t, self.final_t)

class PrincipalCharpoly(Sample):
    """
    characteristic polynomials of principal submatrices
    cauchy interlace theorem
    M is
    5 -4 0
    -4 12 -5
    0 -5 7
    the functions are charpolys of principal submatrices of M.
    f3x3 = t^3 - 24 t^2 + 138 t - 183
    f2x2 = t^2 - 17 t + 44
    f1x1 = t - 5
    """
    def __init__(self):
        t = sympy.abc.t
        p3x3 = sympy.Poly(t**3 - 24*t**2 + 138*t - 183, t)
        p2x2 = sympy.Poly(t**2 - 17*t + 44, t)
        p1x1 = sympy.Poly(t - 5, t)
        roots = [float(r) for r in p3x3.nroots()]
        self.initial_t = min(roots) - 0.05 * (max(roots) - min(roots))
        self.final_t = max(roots) + 0.05 * (max(roots) - min(roots))
        self.polys = (p1x1, p2x2, p3x3)
        self.shape =  interlace.CubicPolyShape(
                self.polys, self.initial_t, self.final_t)
    def get_shape(self):
        return self.shape
    def get_small_3d_sf(self):
        r = self.shape.get_infinity_radius()
        return 0.6 * (1.0 / r)
    def get_large_3d_sf(self):
        return 8.0 * self.get_small_3d_sf()
    def get_axis_radii(self):
        xp_rad = 4.0
        xn_rad = 4.0
        yp_rad = 2.0
        yn_rad = 2.0
        zp_rad = 2.0
        zn_rad = 5.0
        return xp_rad, xn_rad, yp_rad, yn_rad, zp_rad, zn_rad
    def get_superposition_shapes(self):
        return _get_cubic_superposition_shapes(
                self.polys, self.initial_t, self.final_t)

class SchurCharpoly(Sample):
    """
    See sympyutils testing for more info on this specific example.
    characteristic polynomials of schur complement matrices
    """
    def __init__(self):
        t = sympy.abc.t
        p3x3 = sympy.Poly(
                t**3 - 24*(t**2) + 138*t - 183)
        p2x2 = sympy.Poly(
                t**2 - t*sympy.Rational(94, 7) + sympy.Rational(183, 7))
        p1x1 = sympy.Poly(
            t - sympy.Rational(183, 59))
        roots = [float(r) for r in p3x3.nroots()]
        self.initial_t = min(roots) - 0.05 * (max(roots) - min(roots))
        self.final_t = max(roots) + 0.05 * (max(roots) - min(roots))
        self.polys = (p1x1, p2x2, p3x3)
        self.shape =  interlace.CubicPolyShape(
                self.polys, self.initial_t, self.final_t)
    def get_shape(self):
        return self.shape
    def get_small_3d_sf(self):
        r = self.shape.get_infinity_radius()
        return 0.6 * (1.0 / r)
    def get_large_3d_sf(self):
        return 8.0 * self.get_small_3d_sf()
    def get_axis_radii(self):
        xp_rad = 4.0
        xn_rad = 4.0
        yp_rad = 2.0
        yn_rad = 2.0
        zp_rad = 2.0
        zn_rad = 5.0
        return xp_rad, xn_rad, yp_rad, yn_rad, zp_rad, zn_rad
    def get_superposition_shapes(self):
        return _get_cubic_superposition_shapes(
                self.polys, self.initial_t, self.final_t)

class OrthogonalPoly(Sample):
    """
    chebyshev polynomials
    as examples of polynomials that are orthogonal w.r.t. a positive function
    """
    def __init__(self):
        # predefine the shape
        self.initial_t = -1.0
        self.final_t = 1.0
        self.polys = [sympy.chebyshevt_poly(
            i+1, sympy.abc.t, polys=True) for i in range(3)]
        self.shape = interlace.CubicPolyShape(
                self.polys, self.initial_t, self.final_t)
    def get_shape(self):
        return self.shape
    def get_small_3d_sf(self):
        r = self.shape.get_infinity_radius()
        return 0.5 * (1.0 / r)
    def get_large_3d_sf(self):
        return 5.0 * self.get_small_3d_sf()
    def get_axis_radii(self):
        xp_rad = 10.0
        xn_rad = 10.0
        yp_rad = 4.0
        yn_rad = 4.0
        zp_rad = 3.0
        zn_rad = 3.0
        return xp_rad, xn_rad, yp_rad, yn_rad, zp_rad, zn_rad
    def get_superposition_shapes(self):
        return _get_cubic_superposition_shapes(
                self.polys, self.initial_t, self.final_t)

class SturmLiouville(Sample):
    """
    sinusoidal functions
    solutions to -y'' = y with y'=0 at the boundaries
    the functions are cos(n*(t-1)*pi/2)
    as an example of solution of a sturm liouville system
    """
    def __init__(self):
        self.initial_t = -1.0
        self.final_t = 1.0
        self.exprs = [
                sympy.cos(n*(sympy.abc.t-1)*sympy.pi/2) for n in range(1,3+1)]
        self.shape = interlace.DifferentiableShape(
                self.exprs, self.initial_t, self.final_t, 10)
    def get_shape(self):
        return self.shape
    def get_small_3d_sf(self):
        r = self.shape.get_infinity_radius()
        return 0.5 * (1.0 / r)
    def get_large_3d_sf(self):
        return 5.0 * self.get_small_3d_sf()
    def get_axis_radii(self):
        xp_rad = 10.0
        xn_rad = 10.0
        yp_rad = 4.0
        yn_rad = 4.0
        zp_rad = 3.0
        zn_rad = 3.0
        return xp_rad, xn_rad, yp_rad, yn_rad, zp_rad, zn_rad
    def get_superposition_shapes(self):
        return _get_differentiable_superposition_shapes(
                self.exprs, self.initial_t, self.final_t)

class FiniteDifferences(Sample):
    """
    linearly extended eigenvectors of an unweighted path laplacian
    introduce as the finite difference equation solution to
    the sturm liouville system
    use the five segment path
    """
    def __init__(self):
        # define the number of points
        n = 6
        # make the Laplacian matrix
        L = np.zeros((n, n))
        for i in range(n-1):
            L[i, i+1] = -1
            L[i+1, i] = -1
        for i in range(n):
            L[i, i] = 2
        L[0, 0] = 1
        L[-1, -1] = 1
        # define the eigenvector points
        w, vt = scipy.linalg.eigh(L)
        self.x_values = vt.T[1]
        self.y_values = -vt.T[2]
        self.z_values = vt.T[3]
        points = np.array(zip(
            self.x_values, self.y_values, self.z_values))
        edge_lengths = [1.0]*5
        self.shape = interlace.ParametricPiecewiseLinearPathShape(
                points, edge_lengths)
    def get_shape(self):
        return self.shape
    def get_small_3d_sf(self):
        r = self.shape.get_infinity_radius()
        return 0.5 * (1.0 / r)
    def get_large_3d_sf(self):
        return 5.0 * self.get_small_3d_sf()
    def get_axis_radii(self):
        xp_rad = 10.0
        xn_rad = 10.0
        yp_rad = 4.0
        yn_rad = 4.0
        zp_rad = 3.0
        zn_rad = 3.0
        return xp_rad, xn_rad, yp_rad, yn_rad, zp_rad, zn_rad
    def get_superposition_shapes(self):
        shapes = []
        # add the axis shape
        x_axis = interlace.PiecewiseLinearPathShape([
            (0, 0),
            (len(self.x_values)-1, 0)])
        shapes.append(x_axis)
        # add the other piecewise segments
        for values in (self.x_values, self.y_values, self.z_values):
            pairs = list(enumerate(values))
            shape = interlace.PiecewiseLinearPathShape(pairs)
            shapes.append(shape)
        return shapes


class LaplacePath(Sample):
    """
    linearly extended eigenvector of edge-weighted path laplacian
    """
    def __init__(self):
        self.edge_weights = [2.0, 2.0, 3.0, 4.0, 5.0]
        # define the number of points
        n = 6
        # make the Laplacian matrix
        A = np.zeros((n,n))
        for i, weight in enumerate(self.edge_weights):
            A[i+1, i] = weight
            A[i, i+1] = weight
        L = np.diag(A.sum(axis=0)) - A
        # define the eigenvector points
        w, vt = scipy.linalg.eigh(L)
        self.x_values = vt.T[1]
        self.y_values = vt.T[2]
        self.z_values = vt.T[3]
        points = np.array(zip(
            self.x_values, self.y_values, self.z_values))
        edge_lengths = [1.0 / weight for weight in self.edge_weights]
        self.shape = interlace.ParametricPiecewiseLinearPathShape(
                points, edge_lengths)
    def get_shape(self):
        return self.shape
    def get_small_3d_sf(self):
        r = self.shape.get_infinity_radius()
        return 0.5 * (1.0 / r)
    def get_large_3d_sf(self):
        return 5.0 * self.get_small_3d_sf()
    def get_axis_radii(self):
        xp_rad = 10.0
        xn_rad = 10.0
        yp_rad = 4.0
        yn_rad = 4.0
        zp_rad = 3.0
        zn_rad = 3.0
        return xp_rad, xn_rad, yp_rad, yn_rad, zp_rad, zn_rad
    def get_superposition_shapes(self):
        # initialize the times
        t = 0
        times = [t]
        for w in self.edge_weights:
            t += 1.0 / w
            times.append(t)
        # initialize the shapes
        shapes = []
        # add the axis shape
        x_axis = interlace.PiecewiseLinearPathShape([
            (0, 0),
            (times[-1], 0)])
        shapes.append(x_axis)
        # add the other piecewise segments
        for values in (self.x_values, self.y_values, self.z_values):
            pairs = zip(times, values)
            shape = interlace.PiecewiseLinearPathShape(pairs)
            shapes.append(shape)
        return shapes

class LaplaceTree(Sample):
    """
    linearly extended eigenvectors of edge-weighted tree laplacian
    """
    def __init__(self):
        self.tree_string = Newick.daylight_example_tree 
        T, B, N = FtreeIO.newick_to_TBN(self.tree_string)
        self.T = T
        self.B = B
        self.v_to_name = N
        self.shape = self.get_shape()
    def get_small_3d_sf(self):
        r = self.shape.get_infinity_radius()
        return 0.6 * (1.0 / r)
    def get_large_3d_sf(self):
        return 7.0 * self.get_small_3d_sf()
    def get_axis_radii(self):
        xp_rad = 5.0
        xn_rad = 5.0
        yp_rad = 2.0
        yn_rad = 2.0
        zp_rad = 5.0
        zn_rad = 1.5
        return xp_rad, xn_rad, yp_rad, yn_rad, zp_rad, zn_rad
    def get_shape(self):
        # get the full tree laplacian matrix
        vertices = Ftree.T_to_order(self.T)
        L = Ftree.TB_to_L_principal(self.T, self.B, vertices)
        # get the eigendecomposition by increasing eigenvalue
        w, vt = scipy.linalg.eigh(L)
        # get the point valuations of interest
        x_values = vt.T[1]
        y_values = vt.T[2]
        z_values = vt.T[3]
        points = [np.array(xyz) for xyz in zip(x_values, y_values, z_values)]
        # get the vertex to point map
        v_to_point = dict(zip(vertices, points))
        return interlace.PiecewiseLinearTreeShape(self.T, v_to_point)

class SchurTree(Sample):
    """
    harmonically extended eigenvectors
    of schur complement of edge-weighted tree laplacian
    """
    def __init__(self):
        self.tree_string = Newick.daylight_example_tree 
        T, B, N = FtreeIO.newick_to_TBN(self.tree_string)
        self.T = T
        self.B = B
        self.v_to_name = N
        self.shape = self.get_shape()
    def get_small_3d_sf(self):
        r = self.shape.get_infinity_radius()
        return 0.6 * (1.0 / r)
    def get_large_3d_sf(self):
        return 7.0 * self.get_small_3d_sf()
    def get_axis_radii(self):
        xp_rad = 5.0
        xn_rad = 5.0
        yp_rad = 2.0
        yn_rad = 2.0
        zp_rad = 5.0
        zn_rad = 1.5
        return xp_rad, xn_rad, yp_rad, yn_rad, zp_rad, zn_rad
    def get_shape(self):
        # Get the leaf vertices and the internal vertices.
        leaves = Ftree.T_to_leaves(self.T)
        internal = Ftree.T_to_internal_vertices(self.T)
        vertices = leaves + internal
        # Get the harmonic extensions of eigenvectors of schur complement.
        w, v = Ftree.TB_to_harmonic_extension(self.T, self.B, leaves, internal)
        x_values = -v.T[0]
        y_values = -v.T[1]
        z_values = v.T[2]
        points = [np.array(xyz) for xyz in zip(x_values, y_values, z_values)]
        # get the vertex to point map
        v_to_point = dict(zip(vertices, points))
        return interlace.PiecewiseLinearTreeShape(self.T, v_to_point)

