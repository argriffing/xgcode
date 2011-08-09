"""
Draw nine hardcoded logo-sized 3D-embedded 1D shapes.

The command line version will create tikz files.
With basic beamer colors, the default background and foreground are used.
With advanced beamer colors,
the primary and tertiary colors are explicitly pillaged from the theme
and are used to draw the curve and the axes respectively.
Axes are drawn, and the shapes are drawn winding around the axes,
scaled so that the shapes do not go outside the bounding box of the axes.
The orthoplanar intersections are not highlighted
as they are in other renderings.
Because all lines are drawn in the same style,
this simplifies the rendering a little bit.
"""

import argparse
import math

import numpy as np
import sympy

import Form
import FormOut
import tikz
import interlace
import pcurve
import bezier
import sympyutils
import bezintersect

STYLE_AXIS = 'axis-style'
STYLE_AXIS_PATCH = 'axis-patch-style'
STYLE_MAIN = 'main-style'
STYLE_MAIN_PATCH = 'main-patch-style'

def get_form():
    """
    @return: a list of form objects
    """
    # define the form objects
    form_objects = [
            Form.RadioGroup('options', 'color options', [
                Form.RadioItem('black_and_white',
                    'use black and white', True),
                Form.RadioItem('basic_beamer_colors',
                    'use beamer colors'),
                Form.RadioItem('advanced_beamer_colors',
                    'use advanced beamer colors'),
                ]),
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

def get_scene(shape):
    """
    Define all of the bezier paths.
    @return: a list of strokes
    """
    # define the strokes
    strokes = []
    # All of the axis radii are hardcoded.
    # The curve itself is scaled to fit inside the 2d projection of the axes
    xp_rad = xn_rad = 2.0
    yp_rad = yn_rad = 1.0
    zp_rad = zn_rad = 1.0
    strokes.extend([
        bpath_to_stroke(make_half_axis(0, +1, xp_rad), STYLE_AXIS),
        bpath_to_stroke(make_half_axis(0, -1, xn_rad), STYLE_AXIS),
        bpath_to_stroke(make_half_axis(1, +1, yp_rad), STYLE_AXIS),
        bpath_to_stroke(make_half_axis(1, -1, yn_rad), STYLE_AXIS),
        bpath_to_stroke(make_half_axis(2, +1, zp_rad), STYLE_AXIS),
        bpath_to_stroke(make_half_axis(2, -1, zn_rad), STYLE_AXIS)])
    # get the infinity radius of the shape
    r = shape.get_infinity_radius()
    # add the scaled bezier paths of the shape
    for bpath in shape.get_bezier_paths():
        bpath.scale(1.0 / r)
        strokes.append(bpath_to_stroke(bpath, STYLE_MAIN))
    # return the strokes
    return strokes

def get_tikz_pane(shape):
    """
    At this point the tikz styles main-style and axis-style have been defined.
    @param shape: an interlacing.Shape object
    @return: a tikz text string
    """
    min_gridsize = 0.001
    strokes = get_scene(shape)
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
    # get the patches
    patches = []
    for time_list, stroke in zip(time_lists, strokes):
        for patch in stroke.get_patches(time_list):
            if stroke.style == STYLE_MAIN:
                patch.style = STYLE_MAIN_PATCH
            elif stroke.style == STYLE_AXIS:
                patch.style = STYLE_AXIS_PATCH
            patches.append(patch)
    # sort the strokes according to depth order
    depth_stroke_pairs = []
    for stroke in shattered_strokes:
        x, y, z = stroke.evaluate(stroke.characteristic_time)
        depth_stroke_pairs.append((x, stroke))
    # draw the depth sorted strokes and then draw the patches
    arr = []
    ordered_strokes = [s for d, s in sorted(depth_stroke_pairs)] + patches
    for stroke in ordered_strokes:
        # draw a linear curve or a bezier curve
        if len(stroke.bchunks)==1 and stroke.bchunks[0].is_almost_linear():
            p0 = stroke.bchunks[0].p0
            p3 = stroke.bchunks[0].p3
            line = '\\draw[%s] %s -- %s;' % (
                    stroke.style,
                    tikz.point_to_tikz(stroke.bchunks[0].p0[1:]),
                    tikz.point_to_tikz(stroke.bchunks[0].p3[1:]))
            arr.append(line)
        else:
            line = '\\draw[%s]' % stroke.style
            arr.append(line)
            arr.append(get_tikz_bezier(stroke))
    return '\n'.join(arr)


def get_tikz_lines(fs):
    """
    @param fs: user input
    @return: a sequence of tikz lines
    """
    arr = []
    # possibly import some colors from the beamer theme into the tikzpicture
    if fs.advanced_beamer_colors:
        arr.extend([
            '{\\usebeamercolor{palette sidebar tertiary}}',
            '{\\usebeamercolor{palette sidebar primary}}'])
    # determine the foreground and background colors
    if fs.advanced_beamer_colors:
        bg = 'bg'
        fg_main = 'palette sidebar primary.fg'
        fg_axis = 'palette sidebar tertiary.fg'
    elif fs.basic_beamer_colors:
        bg = 'bg'
        fg_main = 'fg'
        fg_axis = 'fg'
    else:
        bg = 'white'
        fg_main = 'black'
        fg_axis = 'black'
    # define the tikzstyles for drawing the curve and the axes
    arr.extend([
        '\\tikzstyle{main-style}=[draw=%s,double=%s,double distance=\\pgflinewidth]' % (bg, fg_main),
        '\\tikzstyle{axis-style}=[draw=%s,double=%s,double distance=\\pgflinewidth]' % (bg, fg_axis),
        '\\tikzstyle{main-patch-style}=[draw=%s]' % fg_main,
        '\\tikzstyle{axis-patch-style}=[draw=%s]' % fg_axis])
    # draw the matrix
    arr.extend([
        '\\matrix{',
        get_tikz_pane(interlace.get_sample_shape_0()),
        '&',
        get_tikz_pane(interlace.get_sample_shape_1()),
        '&',
        get_tikz_pane(interlace.get_sample_shape_2()),
        '\\\\',
        get_tikz_pane(interlace.get_sample_shape_3()),
        '&',
        get_tikz_pane(interlace.get_sample_shape_4()),
        '&',
        get_tikz_pane(interlace.get_sample_shape_5()),
        '\\\\',
        get_tikz_pane(interlace.get_sample_shape_6()),
        '&',
        get_tikz_pane(interlace.get_sample_shape_7()),
        '&',
        get_tikz_pane(interlace.get_sample_shape_8()),
        '\\\\};'])
    return arr

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
    preamble = ''
    document_body = tikz_text
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
