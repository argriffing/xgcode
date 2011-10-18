"""
Draw curve subdivision for two and a half dimensional drawing. [BIT ROTTED]
"""

import math
import random

import numpy as np
import sympy

import Form
import FormOut
import tikz
import interlace
import bezier
import pcurve
import color
import const

def get_form():
    """
    @return: a list of form objects
    """
    # define the form objects
    form_objects = [
            Form.Float('min_gridsize', 'intersection resolution',
                0.1, low_exclusive=0),
            Form.TikzFormat(),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.Tikz()

def rotate_to_view(p):
    """
    Rotate a few degrees around the z axis and then around the new y axis.
    @param p: a 3d point
    @return: a 3d point rotated around the origin
    """
    # use pi/4 for a more standard rotation
    theta = -math.pi / 12
    c = math.cos(theta)
    s = math.sin(theta)
    x0, y0, z0 = p
    # first rotation
    x1 = x0 * c - y0 * s
    y1 = y0 * c + x0 * s
    z1 = z0
    # second rotation
    x2 = x1 * c - z1 * s
    y2 = y1
    z2 = z1 * c + x1 * s
    return np.array([x2, y2, z2])

def get_tikz_lines(fs):
    """
    @param fs: user input
    @return: a sequence of tikz lines
    """
    # define characteristics for all circles
    radius = 1
    nsegments = 5
    # define the first circle
    center = np.array([1.0, 1.0, 1.0])
    axis = 0
    owned_bchunks = bezier.gen_bchunks_ortho_circle(
        center, radius, axis, pcurve.OwnedBezierChunk)
    first_curve = pcurve.BezierPath(owned_bchunks)
    # define the second circle
    center = np.array([1.0, 1.0, 0.0])
    axis = 1
    owned_bchunks = bezier.gen_bchunks_ortho_circle(
        center, radius, axis, pcurve.OwnedBezierChunk)
    second_curve = pcurve.BezierPath(owned_bchunks)
    # rotate every control point in every bchunk in each curve
    for curve in (first_curve, second_curve):
        for b in curve.bchunks:
            b.p0 = rotate_to_view(b.p0)
            b.p1 = rotate_to_view(b.p1)
            b.p2 = rotate_to_view(b.p2)
            b.p3 = rotate_to_view(b.p3)
    # define some new flat curves whose bchunks reference the deep curves
    deep_curves = (first_curve, second_curve)
    flat_curves = []
    for deep_curve in deep_curves:
        flat_bchunks = []
        for deep_b in deep_curve.bchunks:
            flat_b = pcurve.OwnedBezierChunk(
                    deep_b.start_time, deep_b.stop_time,
                    deep_b.p0[1:], deep_b.p1[1:], deep_b.p2[1:], deep_b.p3[1:])
            flat_b.parent_ref = id(deep_curve)
            flat_bchunks.append(flat_b)
        flat_curves.append(pcurve.BezierPath(flat_bchunks))
    # break up the piecewise curves for z-ordering
    child_parent_curve_pairs = list(pcurve.decompose_scene(
            deep_curves, flat_curves, fs.min_gridsize))
    # extract the child curves
    child_curves = zip(*child_parent_curve_pairs)[0]
    # sort the child curves according to depth order
    # TODO chain together the bezier curve into a single drawing command
    depth_curve_pairs = []
    for curve in child_curves:
        x, y, z = rotate_to_view(curve.evaluate(curve.characteristic_time))
        depth_curve_pairs.append((x, curve))
    lines = []
    for x, curve in sorted(depth_curve_pairs):
        colors = ['yellow', 'orange', 'red', 'green', 'blue', 'purple']
        c = random.choice(colors)
        for b in curve.bchunks:
            controls = (b.p0, b.p1, b.p2, b.p3)
            points = [tikz.point_to_tikz(p[1:]) for p in controls]
            args = tuple([c] + points)
            line = '\\draw[draw=white,double=%s,thick] %s .. controls %s and %s .. %s;' % args
            lines.append(line)
            #line = '\\draw %s -- %s -- %s -- %s;' % points
            #lines.append(line)
    return lines

def get_response_content(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: the response
    """
    # get the texts
    tikz_body = '\n'.join(get_tikz_lines(fs))
    tikzpicture = tikz.get_picture(tikz_body, 'auto')
    return tikz.get_response(tikzpicture, fs.tikzformat)
