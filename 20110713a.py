"""
Draw a sequence of superimposed interlacing polynomials.

The input defines a monic cubic polynomial p(t) with distinct real zeros.
"""

import math

import numpy as np

import Form
import FormOut
import tikz
import Ftree
import FtreeIO
import MatrixUtil
import interlace
import pcurve
import const

#TODO finish this

def get_form():
    """
    @return: a list of form objects
    """
    # define the form objects
    form_objects = [
            Form.Float('initial_t', 'initial t', 0.6),
            Form.Float('root_a', 'first root of p(t)', 1.0),
            Form.Float('root_b', 'second root of p(t)', 2.2),
            Form.Float('root_c', 'third root of p(t)', 2.5),
            Form.Float('final_t', 'final t', 3.2),
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

def get_world_segments(root_a, root_b, root_c, initial_t, final_t):
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
    f_x = LineSegment(np.array([-d, 0, 0]), np.array([d, 0, 0]))
    f_y = LineSegment(np.array([0, -d, 0]), np.array([0, d, 0]))
    f_z = LineSegment(np.array([0, 0, -d]), np.array([0, 0, d]))
    x_axis_segs = get_piecewise_curve(f_x, 0, 1, 10, seg_length_min)
    y_axis_segs = get_piecewise_curve(f_y, 0, 1, 10, seg_length_min)
    z_axis_segs = get_piecewise_curve(f_z, 0, 1, 10, seg_length_min)
    segments.extend((p0, p1, STYLE_X) for p0, p1 in x_axis_segs)
    segments.extend((p0, p1, STYLE_Y) for p0, p1 in y_axis_segs)
    segments.extend((p0, p1, STYLE_Z) for p0, p1 in z_axis_segs)
    # add the parametric curve
    f_poly = interlace.InterlacingPoly(
            root_a, root_b, root_c, initial_t, final_t)
    poly_segs = get_piecewise_curve(
            f_poly, initial_t, final_t, 10, seg_length_min)
    segments.extend((p0, p1, STYLE_CURVE) for p0, p1 in poly_segs)
    # add the intersection circles
    radius = 0.1
    x_roots = [f_poly.get_linear_root()]
    y_roots = f_poly.get_quadratic_roots()
    z_roots = f_poly.get_cubic_roots()
    for r in x_roots:
        f = OrthoCircle(f_poly(r), radius, 0)
        segs = get_piecewise_curve(f, 0, 1, 10, seg_length_min)
        segments.extend((p0, p1, STYLE_X) for p0, p1 in segs)
    for r in y_roots:
        f = OrthoCircle(f_poly(r), radius, 1)
        segs = get_piecewise_curve(f, 0, 1, 10, seg_length_min)
        segments.extend((p0, p1, STYLE_Y) for p0, p1 in segs)
    for r in z_roots:
        f = OrthoCircle(f_poly(r), radius, 2)
        segs = get_piecewise_curve(f, 0, 1, 10, seg_length_min)
        segments.extend((p0, p1, STYLE_Z) for p0, p1 in segs)
    # return the segments
    return segments
    

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


def get_tikz_lines(fs):
    """
    @param fs: user input
    @return: a sequence of tikz lines
    """
    segs = get_world_segments(
            fs.root_a, fs.root_b, fs.root_c,
            fs.initial_t, fs.final_t)
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
                STYLE_X: 'red',
                STYLE_Y: 'green',
                STYLE_Z: 'blue',
                STYLE_CURVE: 'black'}[style]
        #line_background = '\\draw[thick,color=white] (%s, %s) -- (%s, %s);' % (y0, z0, y1, z1)
        #line_foreground = '\\draw[color=%s] (%s, %s) -- (%s, %s);' % (color, y0, z0, y1, z1)
        line_double = '\\draw[draw=white,double=%s] (%s, %s) -- (%s, %s);' % (color, y0, z0, y1, z1)
        #lines.append(line_background)
        #lines.append(line_foreground)
        lines.append(line_double)
    # draw dots at the positive endpoints of the axes
    return lines

def get_tikz_text(tikz_body):
    """
    TikZ boilerplate code.
    """
    tikz_header = '\\begin{tikzpicture}[auto]'
    tikz_footer = '\\end{tikzpicture}'
    return '\n'.join([tikz_header, tikz_body, tikz_footer])

def get_latex_text(tikz_text):
    """
    TikZ boilerplate code.
    """
    latex_header = '\n'.join([
        '\\documentclass{article}',
        '\\usepackage{tikz}',
        '\\begin{document}'])
    latex_body = tikz_text
    latex_footer = '\\end{document}'
    return '\n'.join([latex_header, latex_body, latex_footer])

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

