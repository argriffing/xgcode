"""
Create and manipulate one dimensional parametric curves.

The one dimensional parametric curves may live in high dimensional space.
The default embedding space is three dimensional Euclidean space.
Everything in the Bezier section assumes that points are numpy arrays.
"""

from collections import defaultdict
from collections import deque
import unittest
import heapq
import math
import itertools

import numpy as np
from scipy import optimize

import bezier
import iterutils

class BezierPath:
    """
    This curve is created by patching together cubic Bezier curves.
    It may live in a high dimensional space.
    """
    def __init__(self, bchunks):
        """
        @param bchunks: an iterable of BezierChunk objects
        """
        self.bchunks = list(bchunks)
        self.characteristic_time = None
    def get_start_time(self):
        return self.bchunks[0].start_time
    def get_stop_time(self):
        return self.bchunks[-1].stop_time
    def transform(self, f):
        for b in self.bchunks:
            b.transform(f)
    def evaluate_ortho(self, t, axis):
        """
        @param t: time
        @param axis: axis index
        """
        for b in self.bchunks:
            if b.start_time <= t <= b.stop_time:
                return b.eval_global_ortho(t)
    def evaluate(self, t):
        #TODO possibly add a faster function for simultaneous evaluation
        # at multiple times
        for b in self.bchunks:
            if b.start_time <= t <= b.stop_time:
                return b.eval_global(t)
    def get_weak_midpoint_error(self, t_mid, pa, pb):
        p = self.evaluate(t_mid)
        e = np.linalg.norm(p - pa) - np.linalg.norm(pb - p)
        return e*e
    def get_weak_midpoint(self, t_initial, t_final):
        """
        Get a t_mid such that the position is equidistant from the endpoints.
        In other words if c(t) is the curve position at time t,
        then we want a value of t_mid such that
        ||c(t_mid) - c(t_initial)|| == ||c(t_final) - c(t_mid)||.
        Note that this is not necessarily a unique point,
        and it doesn't necessarily have anything to do with arc length.
        @param t_initial: an initial time
        @param t_final: a final time
        @return: t_mid
        """
        args = (self.evaluate(t_initial), self.evaluate(t_final))
        result = scipy.optimize.fminbound(
                self.get_weak_midpoint_error, t_initial, t_final, args)
        return result
    def shatter(self, times):
        """
        Return a collection of BezierPath objects.
        The returned objects should be annotated
        with characteristic times corresponding to intersections.
        @param times: filtered intersection times
        @return: a collection of BezierPath objects
        """
        # handle the edge case of no intersections
        if not times:
            self.characteristic_time = 0.5 * (
                    self.get_start_time() + self.get_stop_time())
            return [self]
        # handle the edge case of a single intersection
        if len(times) == 1:
            self.characteristic_time = times[0]
            return [self]
        # Compute quiescence times.
        # TODO use weak spatially quiescent midpoints
        # instead of naive temporally quiescent midpoints
        quiescence_times = [0.5*(a+b) for a, b in iterutils.pairwise(times)]
        # Construct the bchunks sequences.
        # Use whole bchunks when possible,
        # but at quiescence times we might have to split the bchuncks.
        remaining = deque(self.bchunks)
        groups = []
        g = []
        # repeatedly split the remaining sequence
        for q in quiescence_times:
            while True:
                b = remaining.popleft()
                if b.start_time <= q <= b.stop_time:
                    duration = b.stop_time - b.start_time
                    t_local = (q - b.start_time) / duration
                    alpha, beta = b.split(t_local)
                    g.append(alpha)
                    remaining.appendleft(beta)
                    groups.append(g)
                    g = []
                    break
                else:
                    g.append(b)
        g.extend(remaining)
        groups.append(g)
        # Create a piecewise bezier curve from each group,
        # and give each piecewise curve a characteristic time.
        piecewise_curves = []
        for t, group in zip(times, groups):
            curve = self.__class__(group)
            curve.characteristic_time = t
            piecewise_curves.append(curve)
        return piecewise_curves

def get_bezier_path(fp, fv, t_initial, t_final, nchunks):
    """
    @param fp: a python function from t to position vector
    @param fv: a python function from t to velocity vector
    @param t_initial: initial time
    @param t_final: final time
    @param nchunks: use this many chunks in the piecewise approximation
    @return: a BezierPath
    """
    bchunks = []
    npoints = nchunks + 1
    duration = t_final - t_initial
    incr = duration / nchunks
    times = [t_initial + i*incr for i in range(npoints)]
    for ta, tb in iterutils.pairwise(times):
        b = bezier.create_bchunk_hermite(
                ta, tb, fp(ta), fp(tb), fv(ta), fv(tb))
        bchunks.append(b)
    return BezierPath(bchunks)



#FIXME this function is a linear piecewise holdout that does not use bezier
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

