"""
This is all about Bezier curves.

The curves are assumed to be cubic,
but this code is agnostic to the dimension of the embedding Euclidean space.
The BezierChunk objects of concern in this module are meant to be used
as components of piecewise Bezier curves which should be defined elsewhere.
"""

import math
import itertools

import numpy as np


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

class BezierChunk:
    def __init__(self, start_time, stop_time, p0, p1, p2, p3):
        self.start_time = start_time
        self.stop_time = stop_time
        self.p0 = p0
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
    def eval(self, t_local):
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

