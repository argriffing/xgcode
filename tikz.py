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

def get_picture(tikz_body, options=None):
    """
    @param tikz_body: the text inside a tikzpicture environment
    @param options: a tikzpicture environment options dict 
    @return: the text of a tikzpicture environment
    """
    if options is None:
        options = {}
    options_string = latexutil.options_dict_to_string(options)
    return '\n'.join((
        '\\begin{tikzpicture}' + options_string,
        tikz_body,
        '\\end{tikzpicture}'))

def get_tikz_response(
        packages, preamble, tikz_body, tikzformat, tikzpicture_options=None):
    """
    This is a very simple tikz response.
    For more complicated situations a less generic function
    may be required.
    For example a more specific function may be required
    to add more options to the tikzpicture environment
    or to the imported packages.
    The tikz package is requested automatically.
    @param packages: a collection of requested packages
    @param preamble: color definitions, for example
    @param tikz_body: the text inside a tikzpicture environment
    @param tikzformat: one of four tikz output formats
    @param tikzpicture_options: a tikzpicture environment options dict or None
    @return: a response suitable to return from the get_response interface
    """
    # check the requested format
    assert_tikzformat(tikzformat)
    # define the tikzpicture environment
    if tikzpicture_options is None:
        tikzpicture_options = {'auto' : None}
    tikzpicture = get_picture(tikz_body, tikzpicture_options)
    # immediately return the tikzpicture if requested
    if tikzformat == TIKZFORMAT_TIKZ:
        return tikzpicture
    # delegate to latexutil
    requested_packages = set(packages) | set(['tikz'])
    return latexutil.get_latex_response(
            'standalone', requested_packages, preamble,
            tikzpicture, tikzformat)

