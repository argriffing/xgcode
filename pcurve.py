"""
Create and manipulate one dimensional parametric curves.

The one dimensional parametric curves may live in high dimensional space.
The default embedding space is three dimensional Euclidean space.
"""

import unittest
import heapq
import math

import numpy as np

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

