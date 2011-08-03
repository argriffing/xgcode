"""
Experiment with curve subdivision for two and a half dimensional drawing.
"""

import math
import random

import numpy as np
import sympy

import Form
import FormOut
import tikz
import interlace
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

#TODO remove
def get_world_segments(root_a, root_b, root_c,
        initial_t, final_t, intersection_radius):
    """
    The world consists of
    three axis lines,
    six circles marking the intersections,
    and a single parametric curve.
    @return: a collection of (p0, p1, style) triples
    """
    seg_length_min = 0.1
    segments = []
    # add the axis line segments
    d = 5
    f_x = pcurve.LineSegment(np.array([-d, 0, 0]), np.array([d, 0, 0]))
    f_y = pcurve.LineSegment(np.array([0, -d, 0]), np.array([0, d, 0]))
    f_z = pcurve.LineSegment(np.array([0, 0, -d]), np.array([0, 0, d]))
    x_axis_segs = pcurve.get_piecewise_curve(f_x, 0, 1, 10, seg_length_min)
    y_axis_segs = pcurve.get_piecewise_curve(f_y, 0, 1, 10, seg_length_min)
    z_axis_segs = pcurve.get_piecewise_curve(f_z, 0, 1, 10, seg_length_min)
    segments.extend((p0, p1, STYLE_X) for p0, p1 in x_axis_segs)
    segments.extend((p0, p1, STYLE_Y) for p0, p1 in y_axis_segs)
    segments.extend((p0, p1, STYLE_Z) for p0, p1 in z_axis_segs)
    # add the parametric curve
    roots = (root_a, root_b, root_c)
    polys = interlace.roots_to_differential_polys(roots)
    f_poly = interlace.Multiplex(polys)
    poly_segs = pcurve.get_piecewise_curve(
            f_poly, initial_t, final_t, 10, seg_length_min)
    segments.extend((p0, p1, STYLE_CURVE) for p0, p1 in poly_segs)
    # add the intersection circles
    x_roots_symbolic = sympy.roots(polys[0])
    y_roots_symbolic = sympy.roots(polys[1])
    z_roots_symbolic = sympy.roots(polys[2])
    x_roots = [float(r) for r in x_roots_symbolic]
    y_roots = [float(r) for r in y_roots_symbolic]
    z_roots = [float(r) for r in z_roots_symbolic]
    for r in x_roots:
        f = pcurve.OrthoCircle(f_poly(r), intersection_radius, 0)
        segs = pcurve.get_piecewise_curve(f, 0, 1, 10, seg_length_min)
        segments.extend((p0, p1, STYLE_X) for p0, p1 in segs)
    for r in y_roots:
        f = pcurve.OrthoCircle(f_poly(r), intersection_radius, 1)
        segs = pcurve.get_piecewise_curve(f, 0, 1, 10, seg_length_min)
        segments.extend((p0, p1, STYLE_Y) for p0, p1 in segs)
    for r in z_roots:
        f = pcurve.OrthoCircle(f_poly(r), intersection_radius, 2)
        segs = pcurve.get_piecewise_curve(f, 0, 1, 10, seg_length_min)
        segments.extend((p0, p1, STYLE_Z) for p0, p1 in segs)
    # return the segments
    return segments

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
    first_curve = pcurve.create_bezier_ortho_circle(center, radius, axis)
    # define the second circle
    center = np.array([1.0, 1.0, 0.0])
    axis = 1
    second_curve = pcurve.create_bezier_ortho_circle(center, radius, axis)
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
        curve = pcurve.PiecewiseBezier()
        curve.bchunks = []
        for deep_b in deep_curve.bchunks:
            b = pcurve.BezChunk()
            b.p0 = deep_b.p0[1:]
            b.p1 = deep_b.p1[1:]
            b.p2 = deep_b.p2[1:]
            b.p3 = deep_b.p3[1:]
            b.start_time = deep_b.start_time
            b.stop_time = deep_b.stop_time
            b.parent_ref = id(deep_curve)
            curve.bchunks.append(b)
        flat_curves.append(curve)
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


def get_latex_text(tikz_text):
    """
    TikZ boilerplate code.
    """
    arr = []
    arr.extend([
        '\\documentclass{article}',
        '\\usepackage{tikz}'])
    arr.extend([
        '\\begin{document}',
        tikz_text,
        '\\end{document}'])
    return '\n'.join(arr)

def get_tikz_text(tikz_body):
    """
    TikZ boilerplate code.
    """
    tikz_header = '\\begin{tikzpicture}[auto]'
    tikz_footer = '\\end{tikzpicture}'
    return '\n'.join([tikz_header, tikz_body, tikz_footer])

def get_response_content(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: the response
    """
    # get the texts
    tikz_lines = get_tikz_lines(fs)
    tikz_text = get_tikz_text('\n'.join(tikz_lines))
    latex_text = get_latex_text(tikz_text)
    # decide the output format
    if fs.tikz:
        return tikz_text
    elif fs.tex:
        return latex_text
    elif fs.pdf:
        return tikz.get_pdf_contents(latex_text)
    elif fs.png:
        return tikz.get_png_contents(latex_text)
