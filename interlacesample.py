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
    def get_shape(self):
        return self.shape
    def __init__(self):
        initial_t = 0.9
        root_a = 1.0
        root_b = 2.5
        root_c = 3.5
        final_t = 3.6
        p3 = sympyutils.roots_to_poly((root_a, root_b, root_c))
        p2 = p3.diff()
        p1 = p2.diff()
        polys = (p1, p2, p3)
        self.shape =  interlace.CubicPolyShape(polys, initial_t, final_t)
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

class PrincipalCharpoly(DerivativePoly):
    """
    characteristic polynomials of principal submatrices
    cauchy interlace theorem
    """
    pass

class SchurCharpoly(DerivativePoly):
    """
    characteristic polynomials of schur complement matrices
    """
    pass

class OrthogonalPoly(Sample):
    """
    chebyshev polynomials
    as examples of polynomials that are orthogonal w.r.t. a positive function
    """
    def __init__(self):
        # predefine the shape
        polys = [sympy.chebyshevt_poly(
            i+1, sympy.abc.t, polys=True) for i in range(3)]
        self.shape = interlace.CubicPolyShape(polys, -1.0, 1.0)
    def get_shape(self):
        return self.shape
    def get_small_3d_sf(self):
        r = self.shape.get_infinity_radius()
        return 0.5 * (1.0 / r)
    def get_large_3d_sf(self):
        return 4.0 * self.get_small_3d_sf()
    def get_axis_radii(self):
        xp_rad = 8.0
        xn_rad = 8.0
        yp_rad = 3.0
        yn_rad = 3.0
        zp_rad = 3.0
        zn_rad = 3.0
        return xp_rad, xn_rad, yp_rad, yn_rad, zp_rad, zn_rad

class SturmLiouville(Sample):
    """
    sinusoidal functions
    solutions to -y'' = y with y'=0 at the boundaries
    the functions are cos(n*(t-1)*pi/2)
    as an example of solution of a sturm liouville system
    """
    def __init__(self):
        exprs = [sympy.cos(n*(sympy.abc.t-1)*sympy.pi/2) for n in range(1,3+1)]
        self.shape = interlace.DifferentiableShape(exprs, -1.0, 1.0, 10)
    def get_shape(self):
        return self.shape
    def get_small_3d_sf(self):
        r = self.shape.get_infinity_radius()
        return 0.5 * (1.0 / r)
    def get_large_3d_sf(self):
        return 4.0 * self.get_small_3d_sf()
    def get_axis_radii(self):
        xp_rad = 8.0
        xn_rad = 8.0
        yp_rad = 3.0
        yn_rad = 3.0
        zp_rad = 3.0
        zn_rad = 3.0
        return xp_rad, xn_rad, yp_rad, yn_rad, zp_rad, zn_rad

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
        x_values = vt.T[1]
        y_values = -vt.T[2]
        z_values = vt.T[3]
        points = np.array(zip(*(x_values, y_values, z_values)))
        self.shape = interlace.PiecewiseLinearPathShape(points)
    def get_shape(self):
        return self.shape
    def get_small_3d_sf(self):
        r = self.shape.get_infinity_radius()
        return 0.5 * (1.0 / r)
    def get_large_3d_sf(self):
        return 4.0 * self.get_small_3d_sf()
    def get_axis_radii(self):
        xp_rad = 8.0
        xn_rad = 8.0
        yp_rad = 3.0
        yn_rad = 3.0
        zp_rad = 3.0
        zn_rad = 3.0
        return xp_rad, xn_rad, yp_rad, yn_rad, zp_rad, zn_rad

class LaplacePath(FiniteDifferences):
    """
    linearly extended eigenvector of edge-weighted path laplacian
    """
    pass

class LaplaceTree(FiniteDifferences):
    """
    linearly extended eigenvectors of edge-weighted tree laplacian
    """
    pass

class SchurTree(FiniteDifferences):
    """
    harmonically extended eigenvectors
    of schur complement of edge-weighted tree laplacian
    """
    pass
