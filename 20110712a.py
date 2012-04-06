"""
Draw the parametric plot for a sequence of interlacing polynomials.

The input defines a monic cubic polynomial p(t) with distinct real zeros.
From this cubic polynomial we compute a one dimensional parametric curve
in three dimensional Euclidean space such that f(t) = (p''(t), p'(t), p(t)).
The output is a tikz plot which shows the intersections
between the parametric curve and the planes orthogonal
to the axes of the standard basis.
"""

import math

import numpy as np
import sympy

import Form
import FormOut
import tikz
import interlace
import pcurve
import sympyutils
import color

STYLE_X = 0
STYLE_Y = 1
STYLE_Z = 2
STYLE_CURVE = 3

def get_form():
    """
    @return: a list of form objects
    """
    # define the form objects
    form_objects = [
            Form.Float('initial_t', 'initial t', 0.6),
            Form.Float('root_a', 'first root of p(t)', 1.0),
            Form.Float('root_b', 'second root of p(t)', 2.0),
            Form.Float('root_c', 'third root of p(t)', 2.5),
            Form.Float('final_t', 'final t', 3.2),
            Form.CheckGroup('vis_options', 'visualization options', [
                Form.CheckItem(
                    'fancy_intersect', 'doubly occluded intersections')]),
            Form.Float('radius',
                'curve-plane intersection circle radius', 0.1),
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
    f_poly = interlace.Multiplex((sympyutils.WrappedUniPoly(p) for p in polys))
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
    segs = get_world_segments(
            fs.root_a, fs.root_b, fs.root_c,
            fs.initial_t, fs.final_t, fs.radius)
    # get the rotated and styled line segments
    rotated_and_styled = []
    for p0, p1, style in segs:
        rotated_and_styled.append((
            rotate_to_view(p0),
            rotate_to_view(p1),
            style))
    # do the depth sorting
    quads = []
    for p0, p1, style in rotated_and_styled:
        x0, y0, z0 = p0
        x1, y1, z1 = p1
        depth = min(x0, x1)
        quads.append((depth, tuple(p0), tuple(p1), style))
    # get the tikz lines
    lines = []
    for depth, p0, p1, style in sorted(quads):
        x0, y0, z0 = p0
        x1, y1, z1 = p1
        color = {
                STYLE_X: 'w-blue',
                STYLE_Y: 'w-red',
                STYLE_Z: 'w-olive',
                STYLE_CURVE: 'black'}[style]
        if fs.fancy_intersect:
            line_double = '\\draw[draw=white,double=%s] (%s, %s) -- (%s, %s);' % (color, y0, z0, y1, z1)
            lines.append(line_double)
        else:
            line_foreground = '\\draw[color=%s] (%s, %s) -- (%s, %s);' % (color, y0, z0, y1, z1)
            lines.append(line_foreground)
    # draw dots at the positive endpoints of the axes
    """
    x, y, z = rotate_to_view((3, 0, 0))
    line = '\\draw[fill=red] (%s, %s) circle (0.1);' % (y, z)
    lines.append(line)
    x, y, z = rotate_to_view((0, 3, 0))
    line = '\\draw[fill=green] (%s, %s) circle (0.1);' % (y, z)
    lines.append(line)
    x, y, z = rotate_to_view((0, 0, 3))
    line = '\\draw[fill=blue] (%s, %s) circle (0.1);' % (y, z)
    lines.append(line)
    """
    return lines

def get_response_content(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: the response
    """
    tikz_body = '\n'.join(get_tikz_lines(fs))
    tikzpicture = tikz.get_picture(tikz_body, 'auto')
    return tikz.get_response(
            tikzpicture, fs.tikzformat,
            tikz.get_w_color_package_set(), tikz.get_w_color_preamble())
