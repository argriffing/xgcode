"""
Draw a two and a half dimensional scene with occluding curve intersections.
"""

"""
Draw a scene with occlucing curves using less user input.

The first cubic polynomial root is at t=0.
The third cubic polynomial root is at t=1.
"""

"""
def get_form():
    # define the form objects
    form_objects = [
            Form.Float('root_b',
                'a cubic polynomial root in the interval (0, 1)',
                0.6, low_exclusive=0, high_exclusive=1),
            Form.Float('picture_radius',
                'a tikzpicture radius estimate',
                5, low_exclusive=0),
            Form.TikzFormat(),
            Form.ContentDisposition()]
    return form_objects
"""

import math

import numpy as np
import sympy

import Form
import FormOut
import tikz
import interlace
import pcurve
import bezier
import color
import sympyutils
import bezintersect

STYLE_X = 0
STYLE_Y = 1
STYLE_Z = 2
STYLE_CURVE = 3

g_style_colors = color.wolfram[:3] + ['black']

def get_form():
    """
    @return: a list of form objects
    """
    # define the form objects
    form_objects = [
            Form.Float('initial_t', 'initial t', 0.9),
            Form.Float('root_a', 'first root of p(t)', 1.0),
            Form.Float('root_b', 'second root of p(t)', 2.5),
            Form.Float('root_c', 'third root of p(t)', 3.5),
            Form.Float('final_t', 'final t', 3.6),
            Form.Float('circle_radius',
                'curve-plane intersection circle radius', 0.2, low_exclusive=0),
            Form.Float('x_rad_pos',
                'x+ half-axis radius', 6.0, low_exclusive=0),
            Form.Float('x_rad_neg',
                'x- half-axis radius', 6.0, low_exclusive=0),
            Form.Float('y_rad_pos',
                'y+ half-axis radius', 6.0, low_exclusive=0),
            Form.Float('y_rad_neg',
                'y- half-axis radius', 3.0, low_exclusive=0),
            Form.Float('z_rad_pos',
                'z+ half-axis radius', 3.0, low_exclusive=0),
            Form.Float('z_rad_neg',
                'z- half-axis radius', 3.0, low_exclusive=0),
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

def project_to_2d(point):
    return point[1:]

class Stroke(pcurve.BezierPath):
    def set_style(self, style):
        self.style = style
    def shatter(self, *args, **kwargs):
        pieces = pcurve.BezierPath.shatter(self, *args, **kwargs)
        for p in pieces:
            p.set_style(self.style)
        return pieces

def bpath_to_stroke(bpath, style):
    stroke = Stroke(bpath.bchunks)
    stroke.set_style(style)
    return stroke

def make_half_axis(axis, sign, radius):
    """
    This is a helper function.
    @param axis: in {0, 1, 2}
    @return: bpath
    """
    origin = np.zeros(3)
    target = np.zeros(3)
    target[axis] = sign * radius
    b = bezier.create_bchunk_line_segment(origin, target)
    bpath = pcurve.BezierPath([b])
    b.parent_ref = id(bpath)
    return bpath

def get_scene(root_a, root_b, root_c,
        initial_t, final_t, intersection_radius, half_axis_radii):
    """
    Define all of the bezier paths annotated with styles.
    @return: a list of strokes
    """
    # define the strokes
    strokes = []
    # define the half axis line segments
    xp_rad, xn_rad, yp_rad, yn_rad, zp_rad, zn_rad = half_axis_radii
    strokes.extend([
        bpath_to_stroke(make_half_axis(0, +1, xp_rad), STYLE_X),
        bpath_to_stroke(make_half_axis(0, -1, xn_rad), STYLE_X),
        bpath_to_stroke(make_half_axis(1, +1, yp_rad), STYLE_Y),
        bpath_to_stroke(make_half_axis(1, -1, yn_rad), STYLE_Y),
        bpath_to_stroke(make_half_axis(2, +1, zp_rad), STYLE_Z),
        bpath_to_stroke(make_half_axis(2, -1, zn_rad), STYLE_Z)])
    # define the polynomial curve
    p3 = sympyutils.roots_to_poly((root_a, root_b, root_c))
    p2 = p3.diff()
    p1 = p2.diff()
    shape = interlace.CubicPolyShape((p1, p2, p3), initial_t, final_t)
    strokes.append(bpath_to_stroke(shape.get_bezier_path(), STYLE_CURVE))
    # define the orthocircles at curve-plane intersections
    axes = range(3)
    point_seqs = shape.get_orthoplanar_intersections()
    styles = (STYLE_X, STYLE_Y, STYLE_Z)
    for axis, point_seq, style in zip(axes, point_seqs, styles):
        for center in point_seq:
            bchunks = list(bezier.gen_bchunks_ortho_circle(
                    center, intersection_radius, axis))
            bpath = pcurve.BezierPath(bchunks)
            for b in bchunks:
                b.parent_ref = id(bpath)
            strokes.append(bpath_to_stroke(bpath, style))
    return strokes

def get_tikz_lines(fs):
    """
    @param fs: user input
    @return: a sequence of tikz lines
    """
    lines = []
    min_gridsize = 0.001
    half_axis_radii = (
            fs.x_rad_pos, fs.x_rad_neg,
            fs.y_rad_pos, fs.y_rad_neg,
            fs.z_rad_pos, fs.z_rad_neg)
    strokes = get_scene(
            fs.root_a, fs.root_b, fs.root_c,
            fs.initial_t, fs.final_t,
            fs.circle_radius, half_axis_radii)
    # rotate every control point in every bchunk in each curve
    for stroke in strokes:
        stroke.transform(rotate_to_view)
    # get the intersection times
    time_lists = bezintersect.get_intersection_times(
            strokes, project_to_2d, min_gridsize, 3*min_gridsize)
    # shatter the strokes, tracking the times of interest and the styles
    shattered_strokes = []
    for time_list, stroke in zip(time_lists, strokes):
        shattered_strokes.extend(stroke.shatter(time_list))
    # sort the strokes according to depth order
    depth_stroke_pairs = []
    for stroke in shattered_strokes:
        x, y, z = stroke.evaluate(stroke.characteristic_time)
        depth_stroke_pairs.append((x, stroke))
    # draw the curves
    for x, stroke in sorted(depth_stroke_pairs):
        # draw a linear curve or a bezier curve
        c = g_style_colors[stroke.style]
        if len(stroke.bchunks)==1 and stroke.bchunks[0].is_almost_linear():
            p0 = stroke.bchunks[0].p0
            p3 = stroke.bchunks[0].p3
            line = '\\draw[draw=white,double=%s,thick] %s -- %s;' % (
                    c,
                    tikz.point_to_tikz(stroke.bchunks[0].p0[1:]),
                    tikz.point_to_tikz(stroke.bchunks[0].p3[1:]))
            lines.append(line)
        else:
            line = '\\draw[draw=white,double=%s,thick]' % c
            lines.append(line)
            lines.append(get_tikz_bezier(stroke))
    return lines

def get_tikz_bezier(bpath):
    lines = []
    # draw everything except for the last point of the last chunk
    for b in bpath.bchunks:
        pts = [tikz.point_to_tikz(p[1:]) for p in b.get_points()[:-1]]
        lines.append('%s .. controls %s and %s ..' % tuple(pts))
    # draw the last point of the last chunk
    lines.append('%s;' % tikz.point_to_tikz(bpath.bchunks[-1].p3[1:]))
    return '\n'.join(lines)

def get_latex_text(tikz_text):
    """
    TikZ boilerplate code.
    """
    preamble = '\\usepackage{color}'
    arr = [tikz.define_color(*p) for p in color.wolfram_name_color_pairs]
    arr.append(tikz_text)
    document_body = '\n'.join(arr)
    return tikz.get_latex_text(preamble, document_body)

def get_tikz_text(tikz_body):
    """
    TikZ boilerplate code.
    """
    return '\n'.join([
            '\\begin{tikzpicture}[auto]',
            tikz_body,
            '\\end{tikzpicture}'])

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
