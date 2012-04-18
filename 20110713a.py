"""
Draw a monic cubic polynomial and its derivatives.

The cubic polynomial should have distinct real roots.
"""

import math

import numpy as np
import sympy

import Form
import FormOut
import tikz
import interlace
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

def get_function_tikz_lines(f, t_values, color):
    """
    Note that to use a sympy Poly p you can pass p.eval as f.
    @param p: a univariate function
    @param t_values: a list of values of t
    @param color: something like 'green'
    @return: a list of tikz lines
    """
    arr = ['\\draw[color=%s]' % color]
    i_first = 0
    i_last = len(t_values) - 1
    for i, t in enumerate(t_values):
        if i == i_first:
            prefix = '    '
        else:
            prefix = '    -- '
        if i == i_last:
            suffix = ';'
        else:
            suffix = ''
        line = prefix + str((t, f(t))) + suffix
        arr.append(line)
    return arr

def get_tikz_lines(fs):
    """
    @param fs: user input
    @return: a sequence of tikz lines
    """
    roots = (fs.root_a, fs.root_b, fs.root_c)
    nsegs = 128
    npoints = nsegs + 1
    incr = (fs.final_t - fs.initial_t) / nsegs
    polys = interlace.roots_to_differential_polys(roots)
    t_values = [fs.initial_t + incr*i for i in range(npoints)]
    # get the tikz lines
    lines = []
    lines.extend(get_function_tikz_lines(
        (lambda t: 0), (t_values[0], t_values[-1]), 'black'))
    colors = ('w-blue', 'w-red', 'w-olive')
    for p, color in zip(polys, colors):
        lines.extend(get_function_tikz_lines(
            p.eval, t_values, color))
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
