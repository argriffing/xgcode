"""Help write TikZ things.

This is also for tex scope.
Requires ghostscript and pdflatex.
Default tikz line width is 'thin' which is 0.4pt.
Also pgflinewidth is the current line width.
"""

import tempfile
import subprocess
import os

import iterutils

def get_latex_text(preamble, document_body):
    """
    This may require the tex package called standalone.
    The package may not be packaged for the OS package manager.
    https://help.ubuntu.com/community/LaTeX#Installing%20packages%20manually
    http://www.ctan.org/tex-archive/macros/latex/contrib/standalone
    Besides the .sty I also had to copy the .cls and the .cfg .
    This format is supposed to be compatible with the previewer
    http://www.tlhiv.org/ltxpreview/ .
    """
    return '\n'.join([
        '\\documentclass{standalone}',
        '\\usepackage{tikz}',
        preamble,
        '\\begin{document}',
        '\\thispagestyle{empty}',
        document_body,
        '\\end{document}'])

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
    """
    This is ad hoc and is probably insufficient.
    @param text: unsanitized text
    @return: sanitized text
    """
    arr = []
    d = {
            '\\' : '\\textbackslash{}',
            '_' : '\\textunderscore{}',
            '<' : '\\textless{}',
            '>' : '\\textgreater{}',
            '~' : '\\textasciitilde{}',
            }
    for c in text:
        if c in d:
            arr.append(d[c])
        elif c in '%${}&#':
            arr.append('\\' + c)
        else:
            arr.append(c)
    return ''.join(arr)

def _create_temp_pdf_file(latex_text):
    """
    The returned path name base does not yet have the pdf extension.
    @param latex_text: contents of a LaTeX file
    @return: the base of the path name of a temporary pdf file
    """
    # write a named temporary latex file
    fd, pathname = tempfile.mkstemp(prefix='webtex', dir='/tmp', text=True)
    fout = os.fdopen(fd, 'w+b')
    fout.write(latex_text)
    fout.close()
    # convert the file to a pdf
    args = [
            '/usr/bin/pdflatex',
            '-output-directory', '/tmp',
            '-interaction', 'nonstopmode', pathname]
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()
    p_output = p.stdout.read()
    p_error = p.stderr.read()
    return pathname

def get_png_contents(latex_text):
    """
    This involves using temporary files.
    @param latex_text: contents of a LaTeX file
    @return: contents of a png file
    """
    # create the pdf file
    pathname = _create_temp_pdf_file(latex_text)
    # create the png file
    input_arg = pathname + '.pdf'
    output_arg = '-sOutputFile=%s.png' % pathname
    # sDEVICE used to be pngggray
    args = [
            'gs', '-dSAFER', '-dBATCH', '-dNOPAUSE', '-sDEVICE=png16m',
            output_arg, input_arg]
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()
    p_output = p.stdout.read()
    p_error = p.stderr.read()
    # read the png file
    png_pathname = pathname + '.png'
    try:
        fin = open(png_pathname, 'rb')
        png_contents = fin.read()
        fin.close()
    except IOError, e:
        raise ValueError('failed to create a png file')
    return png_contents

def get_pdf_contents(latex_text):
    """
    This involves using temporary files.
    @param latex_text: contents of a LaTeX file
    @return: contents of a pdf file
    """
    # create the pdf file
    pdf_pathname = _create_temp_pdf_file(latex_text) + '.pdf'
    # read the pdf file
    try:
        fin = open(pdf_pathname, 'rb')
        pdf_contents = fin.read()
        fin.close()
    except IOError, e:
        raise ValueError('failed to create a pdf file')
    return pdf_contents
