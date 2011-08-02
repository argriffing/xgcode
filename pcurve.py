"""
Create and manipulate one dimensional parametric curves.

The one dimensional parametric curves may live in high dimensional space.
The default embedding space is three dimensional Euclidean space.
Everything in the Bezier section assumes that points are numpy arrays.
"""

import unittest
import heapq
import math

import numpy as np

def de_casteljau(p0, p1, p2, p3, t):
    """
    This is a helper function for Bezier evaluation and splitting.
    It is explictly unrolled instead of recursive.
    """
    if not 0 <= t <= 1:
        raise ValueError('t is out of range')
    q01 = (1-t)*p0 + t*p1
    q12 = (1-t)*p1 + t*p2
    q23 = (1-t)*p2 + t*p3
    r012 = (1-t)*q01 + t*q12
    r123 = (1-t)*q12 + t*q23
    s0123 = (1-t)*r012 + t*r123
    return q01, q12, q23, r012, r123, s0123

def bezier_split(p0, p1, p2, p3, t):
    """
    Break the Bezier curve into two new Bezier curves.
    Find ca, cb, cc, cd, ce so that cc is the curve value at t and so that
    the two new bezier curves are (p0, ca, cb, cc) and (cc, cd, ce, p3).
    @return: five new points between p1 and p3
    """
    q01, q12, q23, r012, r123, s0123 = de_casteljau(p0, p1, p2, p3, t)
    return q01, r012, s0123, r123, q23

def bezier_eval(p0, p1, p2, p3, t):
    q01, q12, q23, r012, r123, s0123 = de_casteljau(p0, p1, p2, p3, t)
    return s0123

class BezChunk(object):
    """
    This is a chunk of a Bezier curve for collision testing.
    """
    def __init__(self):
        """
        The refs of original piecewise curves should be hashable and unique.
        I guess that using the object id as a ref should be OK.
        """
        self.parent_ref = None
        self.start_time = None
        self.stop_time = None
        self.p0 = None
        self.p1 = None
        self.p2 = None
        self.p3 = None
    def get_diameter(self):
        return np.linalg.norm(self.get_bb_max() - self.get_bb_min())
    def get_bb_min(self):
        """
        @return: axis aligned bounding box min point
        """
        return np.min([self.p0, self.p1, self.p2, self.p3], axis=0)
    def get_bb_max(self):
        """
        @return: axis aligned bounding box max point
        """
        return np.max([self.p0, self.p1, self.p2, self.p3], axis=0)
    def enumerate_bb_gridpoints(self, gridsize):
        """
        Each yielded gridpoint is a possibly non-unique integer tuple.
        @param gridsize: a positive float
        """
        bb_min = self.get_bb_min() / gridsize
        bb_max = self.get_bb_max() / gridsize
        ranges = []
        for low_float, high_float in zip(bb_min, bb_max):
            low = int(math.floor(low))
            high = int(math.floor(high))
            ranges.append(tuple(range(low, high+1)))
        return itertools.product(*ranges)
    def bisect(self):
        """
        @return: two new BezChunk objects
        """
        t = 0.5
        q01, r012, s0123, r123, q23 = bezier_split(
                self.p0, self.p1, self.p2, self.p3, t)
        # define the first child BezChunk
        a = BezChunk()
        a.parent_ref = self.parent_ref
        a.start_time = self.start_time
        a.stop_time = (1-t)*self.start_time + t*self.stop_time
        a.p0 = self.p0
        a.p1 = q01
        a.p2 = r012
        a.p3 = s0123
        # define the second child BezChunk
        b = BezChunk()
        b.parent_ref = self.parent_ref
        b.start_time = (1-t)*self.start_time + t*self.stop_time
        b.stop_time = self.stop_time
        b.p0 = s0123
        b.p1 = r123
        b.p2 = q23
        b.p3 = self.p3
        # return the new objects
        return a, b


class PiecewiseBezier(object):
    """
    This curve is created by patching together cubic Bezier curves.
    It may live in a high dimensional space.
    """
    def __init__(self):
        pass
    def get_bounding_box(self):
        """
        A Bezier curve falls within the convex hull of its control points.
        Therefore the bounding box of its control points
        is also the bounding box of the curve.
        @return: pmin, pmax
        """

class Bezier(object):
    """
    This is a cubic Bezier curve.
    It may live in a high dimensional space.
    """
    def __init__(self):
        pass
    def get_bounding_box(self):
        """
        @return: pmin, pmax
        """

def get_piecewise_curve(f, t_initial, t_final, npieces_min, seg_length_max):
    """
    Convert a parametric curve into a collection of line segments.
    @param f: returns the (x, y, z) value at time t
    @param t_initial: initial value of t
    @param t_final: final value of t
    @param npieces_min: minimum number of line segments
    @param seg_length_max: maximum line segment length without subdivision
    """
    # define a heap of triples (-length, ta, tb)
    # where length is ||f(tb) - f(ta)||
    q = []
    # initialize the heap
    t_incr = float(t_final - t_initial) / npieces_min
    for i in range(npieces_min):
        ta = t_initial + t_incr * i
        tb = ta + t_incr
        dab = np.linalg.norm(f(tb) - f(ta))
        heapq.heappush(q, (-dab, ta, tb))
    # While segments are longer than the max allowed length,
    # subdivide the segments.
    while -q[0][0] > seg_length_max:
        neg_d, ta, tc = heapq.heappop(q)
        tb = float(ta + tc) / 2
        dab = np.linalg.norm(f(tb) - f(ta))
        dbc = np.linalg.norm(f(tc) - f(tb))
        heapq.heappush(q, (-dab, ta, tb))
        heapq.heappush(q, (-dbc, tb, tc))
    # convert time segments to spatial segments
    return [(f(ta), f(tb)) for neg_d, ta, tb in q]


class OrthoCircle:
    def __init__(self, center, radius, axis):
        """
        @param center: a 3d point
        @param radius: a scalar radius
        @param axis: one of {0, 1, 2}
        """
        self.center = center
        self.radius = radius
        self.axis = axis
    def __call__(self, t):
        """
        @param t: a float in the interval [0, 1]
        @return: a 3d point
        """
        p = np.zeros(3)
        axis_a = (self.axis + 1) % 3
        axis_b = (self.axis + 2) % 3
        theta = 2 * math.pi * t
        p[axis_a] = self.radius * math.cos(theta)
        p[axis_b] = self.radius * math.sin(theta)
        return p + self.center

class LineSegment:
    def __init__(self, initial_point, final_point):
        """
        @param initial_point: the first point of the line segment
        @param final_point: the last point of the line segment
        """
        self.initial_point = initial_point
        self.final_point = final_point
    def __call__(self, t):
        """
        @param t: a float in the interval [0, 1]
        @return: a 3d point
        """
        return self.initial_point * (1-t) + self.final_point * t

