"""
This is all about Bezier curves.

The curves are assumed to be cubic,
but this code is agnostic to the dimension of the embedding Euclidean space.
The BezierChunk objects of concern in this module are meant to be used
as components of piecewise Bezier curves which should be defined elsewhere.
The explicit form of a Bezier curve is
B(t) = (1-t)^3 p0 + 3(1-t)^2 t p1 + 3(1-t) t^2 p2 + t^3 p3
where t is between 0 and 1.
"""

import math
import itertools

import numpy as np
import sympy


class BezierError(ValueError): pass


def de_casteljau(p0, p1, p2, p3, t):
    """
    This is a helper function for Bezier evaluation and splitting.
    It is explictly unrolled instead of recursive.
    """
    if not 0 <= t <= 1:
        raise BezierError('t should be in the interval [0, 1]')
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

def bezier_is_almost_linear(p0, p1, p2, p3, reltol):
    """
    If this returns True then the curve is almost linear.
    If it returns False then it might be linear or might not be linear.
    @param reltol: ratio of displacement to curve length
    """
    linear_distance = np.linalg.norm(p3 - p1)
    if not linear_distance:
        return False
    p1_distance = point_to_line_segment_distance(p1, p0, p3)
    p2_distance = point_to_line_segment_distance(p2, p0, p3)
    return p1_distance + p2_distance < reltol * linear_distance

def point_to_line_segment_distance(x0, x1, x2):
    """
    http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
    @param x0: a test point
    @param x1: an endpoint of the line segment
    @param x2: an endpoint of the line segment
    """
    v01 = x1 - x0
    v12 = x2 - x1
    t = -np.dot(v01, v12) / np.dot(v12, v12)
    if 0 < t < 1:
        # the min distance is to the interior of the line segment
        v = v01 + t*v12
        return np.linalg.norm(v)
    else:
        # the min distance is to an endpoint
        return min(np.linalg.norm(x1 - x0), np.linalg.norm(x2 - x0))

class BezierChunk:
    def __init__(self, start_time, stop_time, p0, p1, p2, p3):
        self.start_time = start_time
        self.stop_time = stop_time
        self.p0 = p0
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
    def is_almost_linear(self, reltol=1e-4):
        return bezier_is_almost_linear(
                self.p0, self.p1, self.p2, self.p3, reltol)
    def global_to_local_time(self, t_global):
        duration = self.stop_time - self.start_time
        return (t_global - self.start_time) / duration
    def eval_global_ortho(self, t_global, axis):
        """
        @param t_global: time in the interval [self.start_time, self.stop_time]
        @param axis: index of the axis of interest
        """
        return self.eval_local_ortho(self.global_to_local_time(t_global))
    def eval_global(self, t_global):
        """
        @param t_global: time in the interval [self.start_time, self.stop_time]
        """
        return self.eval_local(self.global_to_local_time(t_global))
    def eval_local_ortho(self, t_local, axis):
        """
        @param t_local: local time in the interval [0, 1]
        @param axis: index of the axis of interest
        """
        return bezier_eval(
                self.p0[axis], self.p1[axis], self.p2[axis], self.p3[axis],
                t_local)
    def eval_local(self, t_local):
        """
        @param t_local: local time in the interval [0, 1]
        """
        return bezier_eval(self.p0, self.p1, self.p2, self.p3, t_local)
    def get_points(self):
        return self.p0, self.p1, self.p2, self.p3
    def get_diameter(self):
        return np.linalg.norm(self.get_bb_max() - self.get_bb_min())
    def get_bb_min(self):
        """
        @return: axis aligned bounding box min point
        """
        return np.min(self.get_points(), axis=0)
    def get_bb_max(self):
        """
        @return: axis aligned bounding box max point
        """
        return np.max(self.get_points(), axis=0)
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
    def transform(self, f):
        """
        @param f: a transformation function such as a rotation
        """
        self.p0 = f(self.p0)
        self.p1 = f(self.p1)
        self.p2 = f(self.p2)
        self.p3 = f(self.p3)
    def clone(self):
        return self.__class__(
                self.start_time, self.stop_time,
                self.p0, self.p1, self.p2, self.p3)
    def split(self, t):
        """
        @return: two new BezierChunk objects
        """
        q01, r012, s0123, r123, q23 = bezier_split(
                self.p0, self.p1, self.p2, self.p3, t)
        a = self.__class__(
                self.start_time,
                (1-t)*self.start_time + t*self.stop_time,
                self.p0, q01, r012, s0123)
        b = self.__class__(
                (1-t)*self.start_time + t*self.stop_time,
                self.stop_time,
                s0123, r123, q23, self.p3)
        return a, b
    def bisect(self):
        return self.split(0.5)

def create_bchunk_hermite(
        initial_time, final_time,
        initial_point, final_point,
        initial_velocity, final_velocity,
        btype=BezierChunk):
    """
    This function uses the Hermite to Bezier change of basis.
    http://spec.winprog.org/curves/
    @return: a btype object
    """
    duration = final_time - initial_time
    return btype(
            initial_time, final_time,
            initial_point,
            initial_point + (1.0 / 3.0) * (duration * initial_velocity),
            final_point - (1.0 / 3.0) * (duration * final_velocity),
            final_point)

def create_bchunk_line_segment(
        initial_point, final_point, btype=BezierChunk):
    """
    This is a geometric function.
    It assumes that the caller does not care about velocity.
    @return: a btype object
    """
    return btype(
            0.0, 1.0,
            initial_point,
            initial_point * (2.0 / 3.0) + final_point * (1.0 / 3.0),
            initial_point * (1.0 / 3.0) + final_point * (2.0 / 3.0),
            final_point)

def gen_bchunks_ortho_circle(center, radius, axis, btype=BezierChunk):
    """
    This is a geometric function.
    It assumes that the caller does not care about velocity.
    Use a fixed number of segments for the circle.
    Also use a magic number kappa for the control points.
    Yield bchunk objects.
    @param center: a 3d point
    @param radius: a scalar radius
    @param axis: one of {0, 1, 2}
    """
    nsegments = 4
    kappa = 4 * (math.sqrt(2) - 1) / 3
    theta_increment = (math.pi * 2) / nsegments
    # compute the positions and velocities
    axis_a = (axis + 1) % 3
    axis_b = (axis + 2) % 3
    # yield bchunk objects
    for i in range(nsegments):
        # define the angles
        theta_initial = i * theta_increment
        theta_final = (i+1) * theta_increment
        # define the initial point
        p_initial = np.zeros(3)
        p_initial[axis_a] = math.cos(theta_initial)
        p_initial[axis_b] = math.sin(theta_initial)
        # define the final point
        p_final = np.zeros(3)
        p_final[axis_a] = math.cos(theta_final)
        p_final[axis_b] = math.sin(theta_final)
        # define the initial velocity
        v_initial = np.zeros(3)
        v_initial[axis_a] = -math.sin(theta_initial)
        v_initial[axis_b] = math.cos(theta_initial)
        # define the final velocity
        v_final = np.zeros(3)
        v_final[axis_a] = -math.sin(theta_final)
        v_final[axis_b] = math.cos(theta_final)
        # yield the bezier chunk
        yield btype(
                theta_initial, theta_final,
                center + radius * p_initial,
                center + radius * (p_initial + kappa * v_initial),
                center + radius * (p_final - kappa * v_final),
                center + radius * p_final)

