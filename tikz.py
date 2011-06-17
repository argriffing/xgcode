"""Help write TikZ things.

Requires ghostscript and pdflatex.
"""

import tempfile
import subprocess
import os

def ad_hoc_sanitation(text):
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
    args = [
            'gs', '-dSAFER', '-dBATCH', '-dNOPAUSE', '-sDEVICE=pnggray',
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
