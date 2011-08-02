"""
Experiment with Bezier curve intersection.
"""

import math

import numpy as np
import sympy

import Form
import FormOut
import tikz
import color
import pcurve

def get_form():
    """
    @return: a list of form objects
    """
    # define the form objects
    form_objects = [
            Form.Float('min_gridsize', 'intersection resolution',
                0.1, low_exclusive=0),
            Form.TikzFormat(),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.Tikz()

def get_tikz_lines(fs):
    """
    @param fs: user input
    @return: a sequence of tikz lines
    """
    # define a straightforward bezier curve
    p0 = np.array([-1.0, -1.0])
    p1 = np.array([-1.0, 1.0])
    p2 = np.array([1.0, 1.0])
    p3 = np.array([1.0, -1.0])
    # define a messier bezier curve
    q0 = np.array([-2.0, 0.0])
    q1 = np.array([1.0, 1.0])
    q2 = np.array([-1.0, -1.0])
    q3 = np.array([2.0, 0.0])
    # plot the bezier curves using tikz
    arr = []
    points = tuple(tikz.point_to_tikz(p) for p in (p0, p1, p2, p3))
    arr.append('\\draw %s .. controls %s and %s .. %s;' % points)
    points = tuple(tikz.point_to_tikz(p) for p in (q0, q1, q2, q3))
    arr.append('\\draw %s .. controls %s and %s .. %s;' % points)
    # define a BezChunk
    a = pcurve.BezChunk()
    a.p0 = p0
    a.p1 = p1
    a.p2 = p2
    a.p3 = p3
    a.start_time = 0.0
    a.stop_time = 1.0
    a.parent_ref = 10
    # define a BezChunk
    b = pcurve.BezChunk()
    b.p0 = q0
    b.p1 = q1
    b.p2 = q2
    b.p3 = q3
    b.start_time = 0.0
    b.stop_time = 1.0
    b.parent_ref = 11
    # find the intersections
    beziers = pcurve.find_bezier_intersections([a, b], fs.min_gridsize)
    for b in beziers:
        points = tuple(tikz.point_to_tikz(p) for p in (b.p0, b.p1, b.p2, b.p3))
        arr.append('\\draw[red] %s .. controls %s and %s .. %s;' % points)
    # return the lines
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
