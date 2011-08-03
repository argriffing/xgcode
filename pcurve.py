"""
Create and manipulate one dimensional parametric curves.

The one dimensional parametric curves may live in high dimensional space.
The default embedding space is three dimensional Euclidean space.
Everything in the Bezier section assumes that points are numpy arrays.
The two main Bezier classes are BezChunk and PiecewiseBezier.
Objects of these types are called bchunks and curves respectively,
where a curve is basically an aggregation of bchunks.
"""

from collections import defaultdict
from collections import deque
import unittest
import heapq
import math
import itertools

import numpy as np
from scipy import optimize

import iterutils


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
    def gen_bb_gridpoints(self, gridsize):
        """
        Each yielded gridpoint is a possibly non-unique integer tuple.
        @param gridsize: a positive float
        """
        bb_min = self.get_bb_min() / gridsize
        bb_max = self.get_bb_max() / gridsize
        ranges = []
        for low_float, high_float in zip(bb_min, bb_max):
            low = int(math.floor(low_float))
            high = int(math.floor(high_float))
            ranges.append(tuple(range(low, high+1)))
        return itertools.product(*ranges)
    def split(self, t):
        """
        @return: two new BezChunk objects
        """
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
    def bisect(self):
        return self.split(0.5)

def create_bchunk(
        initial_time, final_time,
        initial_point, final_point,
        initial_velocity, final_velocity):
    """
    Create a BezChunk without a parent ref.
    The parent ref may be added later.
    This function uses the Hermite to Bezier change of basis.
    http://spec.winprog.org/curves/
    @return: a BezChunk
    """
    duration = final_time - initial_time
    b = BezChunk()
    b.start_time = initial_time
    b.stop_time = final_time
    b.p0 = initial_point
    b.p1 = initial_point + (1.0 / 3.0) * (duration * initial_velocity)
    b.p2 = final_point - (1.0 / 3.0) * (duration * final_velocity)
    b.p3 = final_point
    return b

def create_bchunk_line_segment(
        initial_time, final_time,
        initial_point, final_point):
    duration = final_time - initial_time
    distance = np.linalg.norm(final_point - initial_point)
    velocity = distance / duration
    b = create_bchunk(
            initial_time, final_time,
            initial_point, final_point,
            velocity, velocity)
    return b

def create_curve_line_segment(
        initial_time, final_time,
        initial_point, final_point):
    """
    Create a PiecewiseBezier consisting of a single BezChunk.
    This curve moves directly from the inital point to the final point
    and the time it takes to get there is equal to the distance.
    @return: a PiecewiseBezier object
    """
    curve = PiecewiseBezier()
    b = create_bchunk_line_segment(
            initial_time, final_time, initial_point, final_point)
    b.parent_ref = id(curve)
    curve.bchunks = [b]
    return curve

def create_curve_ortho_circle(center, radius, axis):
    """
    Use a fixed number of segments for the circle.
    Also use a magic number for the control points.
    @param center: a 3d point
    @param radius: a scalar radius
    @param axis: one of {0, 1, 2}
    @return: a PiecewiseBezier object
    """
    nsegments = 4
    kappa = 4 * (math.sqrt(2) - 1) / 3
    theta_increment = (math.pi * 2) / nsegments
    # compute the positions and velocities
    axis_a = (axis + 1) % 3
    axis_b = (axis + 2) % 3
    # initialize the piecewise bezier curve
    curve = PiecewiseBezier()
    curve.bchunks = []
    for i in range(nsegments):
        # define the angles
        theta_initial = i * theta_increment
        theta_final = (i+1) * theta_increment
        # define the initial point
        p_initial = np.zeros(3)
        p_initial[axis_a] = radius * math.cos(theta_initial)
        p_initial[axis_b] = radius * math.sin(theta_initial)
        p_initial += center
        # define the final point
        p_final = np.zeros(3)
        p_final[axis_a] = radius * math.cos(theta_final)
        p_final[axis_b] = radius * math.sin(theta_final)
        p_final += center
        # define the initial velocity
        v_initial = np.zeros(3)
        v_initial[axis_a] = -radius * math.sin(theta_initial)
        v_initial[axis_b] = radius * math.cos(theta_initial)
        # define the final velocity
        v_final = np.zeros(3)
        v_final[axis_a] = -radius * math.sin(theta_final)
        v_final[axis_b] = radius * math.cos(theta_final)
        # define the bezier chunk
        b = BezChunk()
        b.p0 = p_initial
        b.p1 = p_initial + kappa * v_initial
        b.p2 = p_final - kappa * v_final
        b.p3 = p_final
        b.start_time = theta_initial
        b.stop_time = theta_final
        b.parent_ref = id(curve)
        # add the bezier chunk to the curve
        curve.bchunks.append(b)
    return curve

def create_curve_ortho_circle_naive(center, radius, axis):
    """
    Create a Bezier circle without using the magic number from the internet.
    """
    nsegments = 4
    theta_increment = (math.pi * 2) / nsegments
    # compute the positions and velocities
    axis_a = (axis + 1) % 3
    axis_b = (axis + 2) % 3
    # initialize the piecewise bezier curve
    curve = PiecewiseBezier()
    curve.bchunks = []
    for i in range(nsegments):
        # define the angles
        theta_initial = i * theta_increment
        theta_final = (i+1) * theta_increment
        # define the initial point
        p_initial = np.zeros(3)
        p_initial[axis_a] = radius * math.cos(theta_initial)
        p_initial[axis_b] = radius * math.sin(theta_initial)
        p_initial += center
        # define the final point
        p_final = np.zeros(3)
        p_final[axis_a] = radius * math.cos(theta_final)
        p_final[axis_b] = radius * math.sin(theta_final)
        p_final += center
        # define the initial velocity
        v_initial = np.zeros(3)
        v_initial[axis_a] = -radius * math.sin(theta_initial)
        v_initial[axis_b] = radius * math.cos(theta_initial)
        # define the final velocity
        v_final = np.zeros(3)
        v_final[axis_a] = -radius * math.sin(theta_final)
        v_final[axis_b] = radius * math.cos(theta_final)
        # define the bezier chunk
        b = create_bchunk(
                theta_initial, theta_final,
                p_initial, p_final,
                v_initial, v_final)
        b.parent_ref = id(curve)
        # add the bezier chunk to the curve
        curve.bchunks.append(b)
    return curve

def decompose_scene(deep_curves, flat_curves, min_gridsize):
    """
    The bchunks in the flat curves should reference the deep curves.
    Yield (curve, parent) pairs
    """
    # get the list of all bchunks
    flat_bchunks = []
    for curve in flat_curves:
        flat_bchunks.extend(curve.bchunks)
    # get a smallish set of refined bchunks involved in putative intersections
    intersecting_bchunks = find_bezier_intersections(
            flat_bchunks, min_gridsize)
    # transform the intersecting bchunks into a t set per curve
    curve_id_to_t_set = defaultdict(set)
    for b in intersecting_bchunks:
        curve_id_to_t_set[b.parent_ref].update((b.start_time, b.stop_time))
    # filter the times using the flat curves
    curve_id_to_times = defaultdict(list)
    for flat_curve, deep_curve in zip(flat_curves, deep_curves):
        ref = id(deep_curve)
        t_set = curve_id_to_t_set[ref]
        times = flat_curve.filter_intersection_times(t_set, 3*min_gridsize)
        curve_id_to_times[ref] = times
    # break deep curve into multiple curves
    for curve in deep_curves:
        times = curve_id_to_times.get(id(curve), [])
        child_curves = curve.shatter(times)
        for child in child_curves:
            yield child, curve

def find_bezier_intersections(bchunks, min_gridsize):
    """
    This is essentially a dumb search.
    It looks for collisions of smaller and smaller curve pieces
    on finer and finer grids.
    The only smartness is that if a large curve piece
    has no collisions on a coarse grid,
    then it is not subdivided for consideration in a finer grid search.
    Self intersections of curves are not considered.
    @param bchunks: a collection of BezChunk objects
    @param min_gridsize: a float lower bound resolution
    @return: a collection of refined intersecting bchunks
    """
    # Maintain the invariant that potentially intersecting chunks
    # have a diameter of no more than twice the gridsize.
    gridsize = 0.5 * max(b.get_diameter() for b in bchunks)
    while True:
        # map each grid point to a set of nearby parent curves
        gridmap = defaultdict(set)
        for b in bchunks:
            for gridpoint in b.gen_bb_gridpoints(gridsize):
                gridmap[gridpoint].add(b.parent_ref)
        # Get the set of indices of bchunks
        # whose bounding boxes contain contested grid points.
        index_set = set()
        for i, b in enumerate(bchunks):
            for gridpoint in b.gen_bb_gridpoints(gridsize):
                if len(gridmap[gridpoint]) > 1:
                    index_set.add(i)
                    break
        # Cut the gridsize in half.
        gridsize *= 0.5
        # If the gridsize is below the min
        # then return the bchunks involved in putative intersections.
        if gridsize < min_gridsize:
            return [bchunks[index] for index in index_set]
        # Bisect potentially intersecting chunks
        # until each chunk has a diameter no more than twice the gridsize.
        bchunks_small = []
        bchunks_large = []
        for index in index_set:
            b = bchunks[index]
            if b.get_diameter() <= gridsize*2:
                bchunks_small.append(b)
            else:
                bchunks_large.append(b)
        while bchunks_large:
            b = bchunks_large.pop()
            for child in b.bisect():
                if child.get_diameter() <= gridsize*2:
                    bchunks_small.append(child)
                else:
                    bchunks_large.append(child)
        bchunks = bchunks_small


#TODO possibly add a faster function for simultaneous evaluation
# at multiple times
class PiecewiseBezier(object):
    """
    This curve is created by patching together cubic Bezier curves.
    It may live in a high dimensional space.
    """
    def __init__(self):
        self.bchunks = None
        self.characteristic_time = None
    def get_start_time(self):
        return self.bchunks[0].start_time
    def get_stop_time(self):
        return self.bchunks[-1].stop_time
    def evaluate(self, t):
        for b in self.bchunks:
            if b.start_time <= t <= b.stop_time:
                duration = b.stop_time - b.start_time
                t_local = (t - b.start_time) / duration
                return bezier_eval(b.p0, b.p1, b.p2, b.p3, t_local)
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
    def filter_intersection_times(
            self, raw_intersection_times, min_spatial_gap):
        """
        Collapse intersection clusters.
        @param raw_intersection_times: a collection of intersection times
        @param min_spatial_gap: minimium spatial gap between sequential events
        @return: filtered time sequence
        """
        # first sort the intersection times
        times = sorted(raw_intersection_times)
        # group together times that are indistinguishable
        groups = []
        last_point = None
        g = []
        for t in times:
            point = self.evaluate(t)
            # if we have seen a previous point then check the gap
            if g:
                gap = np.linalg.norm(point - last_point)
                # if the gap is large then start a new group
                if gap >= min_spatial_gap:
                    groups.append(g)
                    g = []
            # append the current time to the current group
            g.append(t)
            # remember the most recent point
            last_point = point
        if g:
            groups.append(g)
        # return the sequence of group midpoints
        return [0.5 * (g[0] + g[-1]) for g in groups]
    def shatter(self, times):
        """
        Return a collection of PiecewiseBezier objects.
        The returned objects should be annotated
        with characteristic times corresponding to intersections.
        @param times: filtered intersection times
        @return: a collection of PiecewiseBezier objects
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
            curve = PiecewiseBezier()
            curve.characteristic_time = t
            curve.bchunks = group
            piecewise_curves.append(curve)
        return piecewise_curves


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

