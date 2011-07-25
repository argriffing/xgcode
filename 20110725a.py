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

def get_form():
    """
    @return: a list of form objects
    """
    # define the form objects
    form_objects = [
            Form.Integer('ncurves', 'draw this many curves',
                3, low=1, high=4),
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
    return -math.cos((math.pi / 2)*(x+1))

def get_tikz_multilines(fs):
    """
    @return: a list of multilines
    """
    ncurves = fs.ncurves
    morph = fs.morph
    # define the number of segments
    nsegs = 128
    npoints = nsegs + 1
    x_initial = -1.0
    x_final = 1.0
    incr = (x_final - x_initial) / float(nsegs)
    # define the colors
    colors = ('wolfram-blue', 'wolfram-red', 'wolfram-olive', 'wolfram-green')
    # define the tikz lines
    arr = []
    # add the black line
    arr.append('\\draw (-1, 0) -- (1, 0);')
    # add the other lines
    for i in range(ncurves):
        mypoly = sympy.chebyshevt_poly(i+1, sympy.abc.x, polys=True)
        points = []
        for j in range(npoints):
            x = x_initial + j*incr
            x_morphed = warp(x) * morph + x * (1 - morph)
            y_morphed = mypoly.eval(x_morphed)
            points.append((x, y_morphed))
        arr.append('\\draw[color=%s]' % colors[i])
        arr.append(tikz.curve_to_tikz(points, 4) + ';')
    return arr

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
    tikz_multilines = get_tikz_multilines(fs)
    tikz_text = get_tikz_text('\n'.join(tikz_multilines))
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

