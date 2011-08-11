"""
Draw nine hardcoded 3D-embedded 1D shapes with color and plane intersections.

The command line version will create tikz files.
Axes are drawn, and the shapes are drawn winding around the axes,
scaled so that the shapes do not go outside the bounding box of the axes.
Intersections with planes orthogonal to the axes are highlighted.
"""

import argparse
import math
import os
import sys

import numpy as np
import sympy

import Form
import FormOut
import tikz
import interlace
import interlacesample
import pcurve
import bezier
import sympyutils
import bezintersect
import color

STYLE_X = 'x-style'
STYLE_X_PATCH = 'x-patch-style'
STYLE_Y = 'y-style'
STYLE_Y_PATCH = 'y-patch-style'
STYLE_Z = 'z-style'
STYLE_Z_PATCH = 'z-patch-style'
STYLE_CURVE = 'curve-style'
STYLE_CURVE_PATCH = 'curve-patch-style'

g_stroke_style_to_patch_style = {
        STYLE_X : STYLE_X_PATCH,
        STYLE_Y : STYLE_Y_PATCH,
        STYLE_Z : STYLE_Z_PATCH,
        STYLE_CURVE : STYLE_CURVE_PATCH}

def get_form():
    """
    @return: a list of form objects
    """
    # define the form objects
    form_objects = [
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

#TODO this might be unused
class Scene:
    """
    This basically has a bunch of (shape, style) pairs.
    """
    def __init__(self, shape_style_pairs):
        self.shape_style_pairs = shape_style_pairs
    def gen_strokes(self):
        for shape, style in self.shape_style_pairs:
            for bpath in shape.get_bezier_paths():
                stroke = Stroke(bpath.bchunks)
                stroke.set_style(style)
                yield stroke


def get_scene(sample):
    """
    Define all of the bezier paths.
    @param sample: an interlacesample.Sample object
    @return: a list of strokes
    """
    # define the strokes
    strokes = []
    # For these full images as opposed to the two-color logos,
    # set the half axis radii per shape.
    # FIXME AD HOC
    xp_rad, xn_rad, yp_rad, yn_rad, zp_rad, zn_rad = sample.get_axis_radii()
    #xp_rad = xn_rad = 3 * (width + height) / 4.0
    #yp_rad = yn_rad = width / 2.0
    #zp_rad = zn_rad = height / 2.0
    # FIXME AD HOC
    # Replace this junk with the right transformations.
    # The right way would involve turning the shapes into
    # bezier paths and then linearly transforming the beziers
    # and then estimating the width and height of the transformation
    # by looking at the bezier bounding boxes and recursively splitting
    # the beziers which are near the x or y max or min, until
    # the x and y max and min are known within some specified error.
    # Then the whole scene can be rescaled to fit within the
    # desired width and height dimensions.
    # The values should be (A^T g) where g is a screen space
    # 2d cardinal direction vector
    # and A is the matrix that rotates and projects the original shapes
    # into screen space.
    """
    bchunks = [b for bpath in shape.get_bezier_paths() for b in bpath.bchunks]
    screen_right = bezier.get_max_dot(bchunks, np.array([-0.3, 0.8, 0]))
    screen_left = bezier.get_max_dot(bchunks, -np.array([-0.3, 0.8, 0]))
    screen_top = bezier.get_max_dot(bchunks, np.array([-0.1, 0.3, 0.8]))
    screen_bottom = bezier.get_max_dot(bchunks, -np.array([-0.1, 0.3, 0.8]))
    scaling_factor = min(
            (width / 2.0) / screen_right,
            (width / 2.0) / screen_left,
            (height / 2.0) / screen_top,
            (height / 2.0) / screen_bottom)
    """
    # define the scaling factor and the shape
    sf = sample.get_large_3d_sf()
    shape = sample.get_shape()
    # add the half axis strokes
    strokes.extend([
        bpath_to_stroke(make_half_axis(0, +1, xp_rad), STYLE_X),
        bpath_to_stroke(make_half_axis(0, -1, xn_rad), STYLE_X),
        bpath_to_stroke(make_half_axis(1, +1, yp_rad), STYLE_Y),
        bpath_to_stroke(make_half_axis(1, -1, yn_rad), STYLE_Y),
        bpath_to_stroke(make_half_axis(2, +1, zp_rad), STYLE_Z),
        bpath_to_stroke(make_half_axis(2, -1, zn_rad), STYLE_Z)])
    # add the scaled bezier paths of the shape
    for bpath in shape.get_bezier_paths():
        bpath.scale(sf)
        strokes.append(bpath_to_stroke(bpath, STYLE_CURVE))
    # define the orthocircles at curve-plane intersections
    intersection_radius = 0.2
    axes = range(3)
    point_seqs = shape.get_orthoplanar_intersections()
    styles = (STYLE_X, STYLE_Y, STYLE_Z)
    for axis, point_seq, style in zip(axes, point_seqs, styles):
        for center in point_seq:
            # NOTE using four-stroke circles
            # NOTE to avoid problems caused by thick back-erased lines
            """
            bchunks = list(bezier.gen_bchunks_ortho_circle(
                    center*sf, intersection_radius, axis))
            bpath = pcurve.BezierPath(bchunks)
            strokes.append(bpath_to_stroke(bpath, style))
            """
            for b in bezier.gen_bchunks_ortho_circle(
                    center*sf, intersection_radius, axis):
                strokes.append(bpath_to_stroke(pcurve.BezierPath([b]), style))
    # return the strokes
    return strokes

def get_tikz_pane(sample):
    """
    At this point the tikz styles main-style and axis-style have been defined.
    @param sample: an interlacesample.Sample object
    @return: a tikz text string
    """
    min_gridsize = 0.001
    strokes = get_scene(sample)
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
    depth_stroke_pairs = []
    for stroke in shattered_strokes:
        x, y, z = stroke.evaluate(stroke.characteristic_time)
        depth_stroke_pairs.append((x, stroke))
    ordered_strokes = [s for d, s in sorted(depth_stroke_pairs)]
    # get the patches, tracking the times of interest and the styles
    patches = []
    for time_list, stroke in zip(time_lists, strokes):
        for patch in stroke.get_patches(time_list):
            patch.style = g_stroke_style_to_patch_style[stroke.style]
            patches.append(patch)
    depth_patch_pairs = []
    for patch in patches:
        x, y, z = patch.evaluate(patch.characteristic_time)
        depth_patch_pairs.append((x, patch))
    ordered_patches = [s for d, s in sorted(depth_patch_pairs)]
    # draw the depth sorted strokes and patches
    arr = []
    #for stroke in ordered_strokes + ordered_patches:
    for stroke in ordered_strokes:
        # draw a linear curve or a bezier curve
        #if len(stroke.bchunks)==1 and stroke.bchunks[0].is_almost_linear():
        if False:
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

def get_tikz_style_definitions():
    """
    return [
            '\\tikzstyle{x-style}=[draw=white,double=w-blue,'
            'double distance=\\pgflinewidth]',
            '\\tikzstyle{y-style}=[draw=white,double=w-red,'
            'double distance=\\pgflinewidth]',
            '\\tikzstyle{z-style}=[draw=white,double=w-olive,'
            'double distance=\\pgflinewidth]',
            '\\tikzstyle{curve-style}=[draw=white,double=black,'
            'double distance=\\pgflinewidth]',
            '\\tikzstyle{x-patch-style}=[draw=w-blue]',
            '\\tikzstyle{y-patch-style}=[draw=w-red]',
            '\\tikzstyle{z-patch-style}=[draw=w-olive]',
            '\\tikzstyle{curve-patch-style}=[draw=black]']
    """
    return [
            '\\tikzstyle{x-style}=[thick,draw=white,double=w-blue,'
            'double distance=\\pgflinewidth]',
            '\\tikzstyle{y-style}=[thick,draw=white,double=w-red,'
            'double distance=\\pgflinewidth]',
            '\\tikzstyle{z-style}=[thick,draw=white,double=w-olive,'
            'double distance=\\pgflinewidth]',
            '\\tikzstyle{curve-style}=[thick,draw=white,double=black,'
            'double distance=\\pgflinewidth]',
            '\\tikzstyle{x-patch-style}=[thick,draw=w-blue]',
            '\\tikzstyle{y-patch-style}=[thick,draw=w-red]',
            '\\tikzstyle{z-patch-style}=[thick,draw=w-olive]',
            '\\tikzstyle{curve-patch-style}=[thick,draw=black]']

def get_tikz_lines():
    """
    @param fs: user input
    @return: a sequence of tikz lines
    """
    arr = []
    # define the tikzstyles for drawing the curve and the axes
    arr.extend(get_tikz_style_definitions())
    # draw the matrix
    samples = interlacesample.get_samples()
    arr.extend([
        '\\matrix{',
        get_tikz_pane(samples[0]),
        '&',
        get_tikz_pane(samples[1]),
        '&',
        get_tikz_pane(samples[2]),
        '\\\\',
        get_tikz_pane(samples[3]),
        '&',
        get_tikz_pane(samples[4]),
        '&',
        get_tikz_pane(samples[5]),
        '\\\\',
        get_tikz_pane(samples[6]),
        '&',
        get_tikz_pane(samples[7]),
        '&',
        get_tikz_pane(samples[8]),
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
    tikz_lines = get_tikz_lines()
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


def main(args):
    for i, sample in enumerate(interlacesample.get_samples()):
        filename = os.path.join(args.outdir, 'sample-%04d.tikz' % i)
        with open(filename, 'w') as fout:
            print 'writing', filename
            arr = []
            # add the remark about the invocation of the generating script
            arr.append('% ' + ' '.join(sys.argv))
            # add the commands to define custom colors
            for name, rgb in color.wolfram_name_color_pairs:
                arr.append(tikz.define_color(name, rgb))
            # add the tikz style definitions
            arr.extend(get_tikz_style_definitions())
            # add the tikz drawing functions
            arr.append(get_tikz_pane(sample))
            # write the file
            print >> fout, '\n'.join(arr)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir',
            default='', help='output directory')
    main(parser.parse_args())

