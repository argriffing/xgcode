"""
Help write TikZ things.

Default tikz line width is 'thin' which is 0.4pt.
Also 'pgflinewidth' is the current line width.
Use 'overlay' if bezier control points are causing trouble;
begin{tikzpicture}[overlay] will stop the bounding box calculation.
"""

import iterutils

class DeprecationError(Exception): pass

def get_latex_text(preamble, document_body, documentclass='standalone'):
    raise InternalError('use latexutil instead')

def point_to_tikz(pt):
    """
    Assume that the format specifier .4f is enough for anybody.
    @param pt: an (x, y) pair
    @return: a tikz string
    """
    return '(' + ', '.join('%.4f' % v for v in pt) + ')'

def curve_to_tikz(points, k=4):
    """
    @param points: (x, y) pairs
    @param k: max number of points per tikz line
    @return: a tikz multiline with up to k points per line
    """
    arr = []
    for chunk in iterutils.ragged_grouper(points, k):
        arr.append(' -- '.join(point_to_tikz(p) for p in chunk))
    return ' --\n'.join(arr)

def define_color(name, rgb):
    """
    @param name: the name of the color
    @param rgb: a triple of integers defining the rgb of the color
    @return: a line of latex code
    """
    r, g, b = rgb
    return '\\definecolor{%s}{RGB}{%s,%s,%s}' % (name, r, g, b)

def sanitize(text):
    raise InternalError('use latexutil instead')

def _create_temp_pdf_file(latex_text):
    raise InternalError('use latexutil instead')

def get_png_contents(latex_text):
    raise InternalError('use latexutil instead')

def get_pdf_contents(latex_text):
    raise InternalError('use latexutil instead')

