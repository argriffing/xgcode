"""
A third visualization of interlacing polynomial sequences.

The cubic polynomial should have distinct real roots.
"""

import math

import numpy as np
import sympy
from sympy.abc import x

import Form
import FormOut
import tikz
import interlace
import iterutils
import const

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
            Form.TikzFormat(),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.Tikz()

def get_segmentation(p, t0, t1):
    """
    A segmentation is a sequence of triples (left, right, sign).
    @param p: a sympy Poly
    @param t0: initial time
    @param t1: final time
    @return: a segmentation
    """
    roots = sorted(float(r) for r in sympy.roots(p))
    points = [t0] + roots + [t1]
    segmentation = []
    for left, right in iterutils.pairwise(points):
        mid = (left + right) / 2
        sign = -1 if p.eval(mid) <= 0 else 1
        seg = (left, right, sign)
        segmentation.append(seg)
    return segmentation

def get_pane_tikz_lines(poly_pair, color_pair, t0, t1):
    """
    This is called maybe three times.
    @param poly_pair: the lower and higher order Poly objects
    @param color_pair: the polynomial colors
    @param t0: initial time
    @param t1: final time
    @return: a list of tikz lines
    """
    lines = []
    pa, pb = poly_pair
    ca, cb = color_pair
    # draw the segmentation corresponding to the lower order polynomial
    for left, right, sign in get_segmentation(pa, t0, t1):
        width = {-1 : '0.5mm', 1 : '2mm'}[sign]
        line = '\\draw[color=%s,line width=%s] (%s, 0) -- (%s, 0);' % (
                ca, width, left, right)
        lines.append(line)
    # draw the cuts corresponding to the higher order polynomial
    for r in sorted(sympy.roots(pb)):
        cut = float(r)
        line = '\\draw[color=%s,line width=0.5mm] (%s, -3mm) -- (%s, 3mm);' % (
                cb, cut, cut)
        lines.append(line)
    return lines

def get_tikz_lines(fs):
    """
    @param fs: user input
    @return: a sequence of tikz lines
    """
    # construct a 'matrix' of three vertically stacked panes
    colors = ('black', 'red', 'green', 'blue')
    p0 = sympy.Poly(1, sympy.abc.x)
    roots = (fs.root_a, fs.root_b, fs.root_c)
    polys = [p0] + interlace.roots_to_differential_polys(roots)
    # get the tikz lines
    lines = []
    lines.append('\\matrix {')
    lines.extend(get_pane_tikz_lines(
        polys[0:2], colors[0:2], fs.initial_t, fs.final_t))
    lines.append('\\\\')
    lines.extend(get_pane_tikz_lines(
        polys[1:3], colors[1:3], fs.initial_t, fs.final_t))
    lines.append('\\\\')
    lines.extend(get_pane_tikz_lines(
        polys[2:4], colors[2:4], fs.initial_t, fs.final_t))
    lines.append('\\\\')
    lines.append('};')
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
    Add the matrix library.
    """
    latex_header = '\n'.join([
        '\\documentclass{article}',
        '\\usepackage{tikz}',
        '\\usetikzlibrary{matrix}',
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

