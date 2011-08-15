"""
Draw nine hardcoded interlacing shape superpositions.

The command line version will create tikz files.
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
import typeutils

STYLE_X = 'x-style'
STYLE_Y = 'y-style'
STYLE_Z = 'z-style'
STYLE_CURVE = 'curve-style'


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

def get_tikz_pane(sample, width=6, height=6):
    shapes = None
    try:
        shapes = sample.get_superposition_shapes()
    except AttributeError, e:
        pass
    if shapes:
        return interlace.tikz_shape_superposition(shapes, width, height)
    else:
        return _get_tikz_pane_tree(sample)

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

def get_tikz_style_definitions():
    return [
            '\\tikzstyle{x-style}=[thick,draw=white,double=w-blue,'
            'double distance=\\pgflinewidth]',
            '\\tikzstyle{y-style}=[thick,draw=white,double=w-red,'
            'double distance=\\pgflinewidth]',
            '\\tikzstyle{z-style}=[thick,draw=white,double=w-olive,'
            'double distance=\\pgflinewidth]',
            '\\tikzstyle{curve-style}=[thick,draw=white,double=black,'
            'double distance=\\pgflinewidth]']

def get_tikz_bezier(bpath):
    lines = []
    # draw everything except for the last point of the last chunk
    for b in bpath.bchunks:
        pts = [tikz.point_to_tikz(p[1:]) for p in b.get_points()[:-1]]
        lines.append('%s .. controls %s and %s ..' % tuple(pts))
    # draw the last point of the last chunk
    lines.append('%s;' % tikz.point_to_tikz(bpath.bchunks[-1].p3[1:]))
    return '\n'.join(lines)

def _get_scene(tree_sample):
    """
    Define all of the bezier paths.
    @param tree_sample: an interlacesample.Sample object
    @return: a list of strokes
    """
    # define the strokes
    strokes = []
    # define the scaling factor and the shape
    sf = 1.0
    styles = (STYLE_CURVE, STYLE_X, STYLE_Y, STYLE_Z)
    layout_width = 8.0
    layout_height = 8.0
    value_height = 8.0
    shapes = tree_sample.get_tree_superposition_shapes(
            layout_width, layout_height, value_height)
    for shape, style in zip(shapes, styles):
        # add the scaled bezier paths of the shape
        for bpath in shape.get_bezier_paths():
            bpath.scale(sf)
            strokes.append(bpath_to_stroke(bpath, style))
    # return the strokes
    return strokes

def _get_tikz_pane_tree(tree_sample):
    """
    @return: a tikz text string
    """
    min_gridsize = 0.001
    strokes = _get_scene(tree_sample)
    # rotate every control point in every bchunk in each curve
    for stroke in strokes:
        stroke.transform(interlacesample.rotate_to_view)
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
    # draw the depth sorted strokes and patches
    arr = []
    for stroke in ordered_strokes:
        line = '\\draw[%s]' % stroke.style
        arr.append(line)
        arr.append(get_tikz_bezier(stroke))
    return '\n'.join(arr)

def get_tikz_lines():
    """
    @param fs: user input
    @return: a sequence of tikz lines
    """
    arr = []
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

def get_latex_text(tikz_text):
    """
    TikZ boilerplate code.
    """
    preamble = '\\usepackage{color}'
    # add color definitions
    arr = [tikz.define_color(*p) for p in color.wolfram_name_color_pairs]
    # add style definitions for the trees
    arr.extend(get_tikz_style_definitions())
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
            # add style definitions for the trees
            arr.extend(get_tikz_style_definitions())
            # add the tikz drawing functions
            arr.append(get_tikz_pane(sample, args.width, args.height))
            # write the file
            print >> fout, '\n'.join(arr)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--width',
            default=6, type=typeutils.positive_float,
            help='max width in default tikz units')
    parser.add_argument('--height',
            default=6, type=typeutils.positive_float,
            help='max height in default tikz units')
    parser.add_argument('--outdir',
            default='', help='output directory')
    main(parser.parse_args())

