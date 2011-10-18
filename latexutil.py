"""
Utilities for LaTeX.

This will probably work only for texlive.
Also it requires ghostscript to make a png file.
"""

import tempfile
import subprocess
import os

import iterutils

class CheckPackageError(Exception): pass

def check_packages(package_names):
    """
    @param package_names: a collection of package names
    @return: the installed subset of package names
    """
    # get the file names corresponding to the packages
    filename_args = [name + '.sty' for name in package_names]
    # run the kpsewhich program and capture its output
    args = ['kpsewhich'] + filename_args
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()
    p_output = p.stdout.read()
    p_error = p.stderr.read()
    # get the subset of installed package names
    requested_name_set = set(package_names)
    installed_name_set = set()
    for raw_line in p_output.splitlines():
        line = raw_line.strip()
        if line.endswith('.sty'):
            name, ext = os.path.splitext(os.path.basename(line))
            if name in requested_name_set:
                installed_name_set.add(name)
            else:
                raise CheckPackageError(line)
    # return the set of installed package names
    return installed_name_set

def get_latex_text(preamble, document_body, documentclass='standalone'):
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
        '\\documentclass{%s}' % documentclass,
        '\\usepackage{tikz}',
        preamble,
        '\\begin{document}',
        '\\thispagestyle{empty}',
        document_body,
        '\\end{document}'])

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
            '-output-directory=/tmp',
            '-interaction=nonstopmode',
            '-halt-on-error',
            pathname]
    # 
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
    png_pathname = pathname + '.png'
    input_arg = pathname + '.pdf'
    output_arg = '-sOutputFile=' + png_pathname
    # sDEVICE used to be pngggray
    args = [
            'gs', '-dSAFER', '-dBATCH', '-dNOPAUSE', '-sDEVICE=png16m',
            output_arg, input_arg]
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()
    p_output = p.stdout.read()
    p_error = p.stderr.read()
    # read the png file
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

