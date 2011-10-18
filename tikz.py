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

g_tikzformats = set([
    TIKZFORMAT_TIKZ,
    TIKZFORMAT_TEX,
    TIKZFORMAT_PDF,
    TIKZFORMAT_PNG])

class DeprecationError(Exception): pass
class LatexPackageError(Exception): pass


def get_latex_text(preamble, document_body, documentclass='standalone'):
    raise DeprecationError('use latexutil instead')

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

def get_tikz_response(packages, preamble, tikz_body, tikzformat):
    """
    This is a very generic tikz response.
    For more complicated situations a less generic function
    may be required.
    For example a more specific function may be required
    to add more options to the tikzpicture environment
    or to the imported packages.
    The tikz package is requested automatically.
    @param packages: a collection of requested packages
    @param preamble: color definitions, for example
    @param tikz_body: tikzpicture contents
    @param tikzformat: one of four tikz output formats
    @return: a response suitable to return from the get_response interface
    """
    # check the requested format
    if tikzformat not in g_tikzformats:
        msg = 'invalid requested format: ' + tikzformat
        raise ValueError(msg)
    # define some sets of requested packages
    user_requested_set = set(packages) | set(['tikz'])
    requested_set = user_requested_set | set(['standalone'])
    # define the subset of available packages
    installed_set = set(latexutil.check_packages(requested_set))
    # If no essential packages are missing
    # then define the documentclass and usepackage lines.
    missing_list = sorted(user_requested_set - installed_set)
    if not missing_list:
        if 'standalone' in installed_set:
            documentclass = 'standalone'
        else:
            documentclass = 'article'
        installed_list = sorted(installed_set)
        package_lines = ['\\usepackage{%s}' % s for s in installed_list]
    # If essential packages are missing and a compiled response
    # has been requested then this is a problem.
    # If essential packages are missing but a compiled response
    # has not been requested, then show the most optimistic settings.
    if missing_list:
        if tikzformat in (TIKZFORMAT_PDF, TIKZFORMAT_PNG):
            msg = 'missing LaTeX packages: ' + ' '.join(missing_list)
            raise LatexPackageError(msg)
        else:
            documentclass = 'standalone'
            requested_list = sorted(requested_set)
            package_lines = ['\\usepackage{%s}' % s for s in requested_list]
    # define the latex body
    latex_text = '\n'.join((
        '\\documentclass{%s}' % documentclass,
        '\n'.join(package_lines),
        preamble,
        '\\begin{document}',
        '\\begin{tikzpicture}[auto]',
        tikz_body,
        '\\end{tikzpicture}',
        '\\end{document}'))
    # respond using the requested format
    if tikzformat == TIKZFORMAT_TIKZ:
        return tikz_body
    elif tikzformat == TIKZFORMAT_TEX:
        return latex_text
    elif tikzformat == TIKZFORMAT_PDF:
        return latexutil.get_pdf_contents(latex_text)
    elif tikzformat == TIKZFORMAT_PNG:
        return latexutil.get_png_contents(latex_text)

def _create_temp_pdf_file(latex_text):
    raise DeprecationError('use latexutil instead')

def get_png_contents(latex_text):
    raise DeprecationError('use latexutil instead')

def get_pdf_contents(latex_text):
    raise DeprecationError('use latexutil instead')

