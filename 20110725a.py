"""
Draw Chebyshev polynomials stretched into sinusoidal functions.
"""

import math

import sympy
import sympy.abc

import Form
import FormOut
import tikz
import color
import interlace

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

def get_tikzpicture_body(fs):
    """
    Try to use sympy for computation and bezier for drawing.
    """
    ncurves = fs.ncurves
    morph = fs.morph
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
                (x_expr, y_expr), -1.0, 1.0, fs.nsegs)
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

def get_latex_text(tikz_text):
    """
    TikZ boilerplate code.
    """
    arr = []
    arr.extend([
        '\\documentclass{article}',
        '\\usepackage{tikz}',
        '\\usepackage{color}'])
    arr.extend(
        tikz.define_color(*pair) for pair in color.wolfram_name_color_pairs)
    arr.extend([
        '\\begin{document}',
        tikz_text,
        '\\end{document}'])
    return '\n'.join(arr)

def get_tikzpicture(tikz_body):
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
    tikz_text = get_tikzpicture(get_tikzpicture_body(fs))
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

