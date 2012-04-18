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
import color

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
            Form.TikzFormat()]
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
    colors = ('black', 'w-blue', 'w-red', 'w-olive')
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
