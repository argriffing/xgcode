"""
Help write TikZ things.

Default tikz line width is 'thin' which is 0.4pt.
Also 'pgflinewidth' is the current line width.
Use 'overlay' if bezier control points are causing trouble;
begin{tikzpicture}[overlay] will stop the bounding box calculation.
The nice way to make standalone pdfs which contain only a single
tikz picture is to use the tex package called 'standalone'.
The package may not be packaged for the OS package manager.
https://help.ubuntu.com/community/LaTeX#Installing%20packages%20manually
http://www.ctan.org/tex-archive/macros/latex/contrib/standalone
Besides the .sty I also had to copy the .cls and the .cfg .
This format is supposed to be compatible with the previewer
http://www.tlhiv.org/ltxpreview/ .
"""

import iterutils
import latexutil
import color

TIKZFORMAT_TIKZ = 'tikz'
TIKZFORMAT_TEX = 'tex'
TIKZFORMAT_PDF = 'pdf'
TIKZFORMAT_PNG = 'png'

g_tikzformats = set((
    TIKZFORMAT_TIKZ,
    TIKZFORMAT_TEX,
    TIKZFORMAT_PDF,
    TIKZFORMAT_PNG))


# This is an example from texample
# \usepackage{tikz}
# \usepackage{verbatim}
g_stub = r"""
% Draw axes
\draw [<->,thick] (0,2) node (yaxis) [above] {$y$}
|- (3,0) node (xaxis) [right] {$x$};
% Draw two intersecting lines
\draw (0,0) coordinate (a_1) -- (2,1.8) coordinate (a_2);
\draw (0,1.5) coordinate (b_1) -- (2.5,0) coordinate (b_2);
% Calculate the intersection of the lines a_1 -- a_2 and b_1 -- b_2
% and store the coordinate in c.
\coordinate (c) at (intersection of a_1--a_2 and b_1--b_2);
% Draw lines indicating intersection with y and x axis. Here we use
% the perpendicular coordinate system
\draw[dashed] (yaxis |- c) node[left] {$y'$}
-| (xaxis -| c) node[below] {$x'$};
% Draw a dot to indicate intersection point
\fill[red] (c) circle (2pt);
""".strip()

def assert_tikzformat(tikzformat):
    if tikzformat not in g_tikzformats:
        msg = 'invalid requested format: ' + tikzformat
        raise ValueError(msg)

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

def get_w_color_preamble():
    arr = [define_color(*p) for p in color.wolfram_name_color_pairs]
    color_preamble = '\n'.join(arr)
    return color_preamble

def get_w_color_package_set():
    return set(['color'])

def get_picture(tikz_body, *args, **kwargs):
    """
    @param tikz_body: the text inside a tikzpicture environment
    @param args: single tikzpicture options
    @param kwargs: paired tikzpicture options
    @return: the text of a tikzpicture environment
    """
    redundant_pair = ('scale', 1)
    kwargs = dict(p for p in kwargs.items() if p != redundant_pair)
    option_string = latexutil.options_to_string(*args, **kwargs)
    return '\n'.join((
        '\\begin{tikzpicture}' + option_string,
        tikz_body,
        '\\end{tikzpicture}'))

def get_response(tikzpicture, tikzformat, packages=(), preamble=''):
    """
    @param tikzpicture: a complete tikzpicture environment
    @param tikzformat: one of four tikz output formats
    @param packages: a collection of requested packages
    @param preamble: color definitions, for example
    @return: a response suitable to return from the get_response interface
    """
    # check the requested format
    assert_tikzformat(tikzformat)
    # immediately return the tikzpicture if requested
    if tikzformat == TIKZFORMAT_TIKZ:
        return tikzpicture
    # delegate to latexutil
    requested_packages = set(packages) | set(['tikz'])
    return latexutil.get_response(
            'standalone', tikzpicture, tikzformat,
            requested_packages, preamble)

def get_figure_response(
        tikzpicture, tikzformat, figure_caption, figure_label,
        packages=(), preamble=''):
    """
    @param tikzpicture: a complete tikzpicture environment
    @param tikzformat: one of four tikz output formats
    @param figure_caption: figure caption
    @param figure_label: figure label
    @param packages: a collection of requested packages
    @param preamble: color definitions, for example
    @return: a response suitable to return from the get_response interface
    """
    # check the requested format
    assert_tikzformat(tikzformat)
    # immediately return the tikzpicture if requested
    if tikzformat == TIKZFORMAT_TIKZ:
        return tikzpicture
    # delegate to latexutil
    requested_packages = set(packages) | set(['tikz'])
    return latexutil.get_centered_figure_response(
            tikzpicture, tikzformat, figure_caption, figure_label,
            requested_packages, preamble)

