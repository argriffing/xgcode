"""
Draw Chebyshev polynomials stretched into sinusoidal functions.

Use the command line to spam some morphs into separate tikz files
to be included in a larger tex document.
"""

import argparse
import math
import os
import sys

import sympy
import sympy.abc

import Form
import FormOut
import tikz
import color
import interlace
import typeutils

def get_form():
    """
    @return: a list of form objects
    """
    # define the form objects
    form_objects = [
            Form.Integer('ncurves', 'draw this many curves',
                3, low=1, high=4),
            Form.Integer('nsegs', 'use this many beziers per curve',
                10, low=1, high=100),
            Form.Float('morph', 'morph progress in [0,1]',
                0.5, low_inclusive=0, high_inclusive=1),
            Form.TikzFormat(),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.Tikz()

def warp(x):
    """
    @param x: a scalar in [-1, 1]
    @return: another scalar in [-1, 1]
    """
    return math.sin((math.pi / 2)*x)

def inv_warp(x):
    """
    @param x: a scalar in [-1, 1]
    @return: another scalar in [-1, 1]
    """
    return (2 / math.pi) * math.asin(x)

def get_tikzpicture_body(ncurves, nsegs, morph):
    """
    Try to use sympy for computation and bezier for drawing.
    """
    # define the sympy expression for the warp
    sympy_t = sympy.abc.t
    sympy_t_warped = sympy.sin((sympy.pi / 2) * sympy_t)
    sympy_t_morphed = morph*sympy_t_warped + (1-morph)*sympy_t
    # define the shapes
    shapes = []
    x_axis = interlace.PiecewiseLinearPathShape([(-1, 0), (1, 0)])
    shapes.append(x_axis)
    for i in range(ncurves):
        x_expr = sympy_t
        y_expr = sympy.chebyshevt_poly(i+1, sympy_t_morphed)
        #cheby_expr = sympy.chebyshevt_poly(i+1, sympy_t)
        #y_expr = cheby_expr.subs(sympy_t, sympy_t_morphed)
        shape = interlace.DifferentiableShape(
                (x_expr, y_expr), -1.0, 1.0, nsegs)
        shapes.append(shape)
    width = 6
    height = 6
    return interlace.tikz_shape_superposition(shapes, width, height)

def get_tikzpicture_body_gamma(fs):
    """
    This was the latest pre-bezier version.
    """
    ncurves = fs.ncurves
    morph = fs.morph
    # define the number of segments
    nsegs = 128
    npoints = nsegs + 1
    t_initial = -1.0
    t_final = 1.0
    incr = (t_final - t_initial) / float(nsegs)
    # define the sequence of t values
    t_seq = [t_initial + t*incr for t in range(npoints)]
    # define the sequences of y values
    y_seqs = []
    for i in range(ncurves):
        mypoly = sympy.chebyshevt_poly(i+1, sympy.abc.x, polys=True)
        y_seq = [mypoly.eval(warp(t)*morph + t*(1-morph)) for t in t_seq]
        y_seqs.append(y_seq)
    width = 8
    height = 4
    return interlace.tikz_superposition(t_seq, y_seqs, width, height)

def get_tikzpicture_body_beta(fs):
    """
    This alternative visualization literally warps the x axis.
    The other one makes a smoother picture by
    computing a different function for y that is equivalent to this warping.
    """
    ncurves = fs.ncurves
    morph = fs.morph
    # define the number of segments
    nsegs = 128
    npoints = nsegs + 1
    t_initial = -1.0
    t_final = 1.0
    incr = (t_final - t_initial) / float(nsegs)
    # define the sequence of t values
    t_seq = [t_initial + t*incr for t in range(npoints)]
    t_seq_inv_warped = [
            inv_warp(t)*morph + t*(1-morph) for t in t_seq]
    # define the sequences of y values
    y_seqs = []
    for i in range(ncurves):
        mypoly = sympy.chebyshevt_poly(i+1, sympy.abc.x, polys=True)
        y_seqs.append([mypoly.eval(t) for t in t_seq])
    width = 8
    height = 4
    text = interlace.tikz_superposition(
            t_seq_inv_warped, y_seqs, width, height)
    return text

def get_response_content(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: the response
    """
    tikz_body = get_tikzpicture_body(fs.ncurves, fs.nsegs, fs.morph)
    tikzpicture = tikz.get_picture(tikz_body, 'auto')
    return tikz.get_response(
            tikzpicture, fs.tikzformat,
            tikz.get_w_color_package_set(), tikz.get_w_color_preamble())

def main(args):
    for i in range(args.nframes):
        morph = i / float(args.nframes - 1)
        filename = os.path.join(args.outdir, 'frame%04d.tikz' % i)
        with open(filename, 'w') as fout:
            print 'writing', filename
            arr = []
            # add the remark about the invocation of the generating script
            arr.append('% ' + ' '.join(sys.argv))
            # add the color definitions
            for name, rgb in color.wolfram_name_color_pairs:
                arr.append(tikz.define_color(name, rgb))
            # add the tikz
            arr.append(get_tikzpicture_body(args.ncurves, args.nsegs, morph))
            print >> fout, '\n'.join(arr)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ncurves', type=typeutils.positive_integer,
            default=3, help='plot this many curves on top of the x axis')
    parser.add_argument('--nsegs', type=typeutils.positive_integer,
            default=10, help='use a piecewise bezier with this many segments')
    parser.add_argument('--nframes', type=typeutils.positive_integer,
            default=3, help='make this many files like frame????.tikz')
    parser.add_argument('--outdir',
            default='', help='output directory')
    main(parser.parse_args())

