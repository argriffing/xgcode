"""
Utilities for LaTeX.

This will probably work only for texlive.
Also it requires ghostscript or imagemagick to make a png file.
It seems like the png files made by ghostscript are not so great.
"""

import unittest
import tempfile
import subprocess
import os


LATEXFORMAT_TEX = 'tex'
LATEXFORMAT_PDF = 'pdf'
LATEXFORMAT_PNG = 'png'

g_latexformats = {LATEXFORMAT_TEX, LATEXFORMAT_PDF, LATEXFORMAT_PNG}


class CheckPackageError(Exception): pass
class LatexPackageError(Exception): pass


def assert_latexformat(latexformat):
    if latexformat not in g_latexformats:
        msg = 'invalid requested format: ' + latexformat
        raise ValueError(msg)

def options_to_string(*args, **kwargs):
    """
    Everything is stringified.
    @param args: single values
    @param kwargs: argument pairs
    @return: a string
    """
    arr = []
    for v in sorted(args):
        arr.append(v)
    for k, v in sorted(kwargs.items()):
        arr.append(str(k) + '=' + str(v))
    if arr:
        return '[' + ','.join(arr) + ']'
    else:
        return ''

def _check_installed_files(filenames):
    """
    @param filenames: a collection of requested filenames
    @return: the subset of installed filenames
    """
    requested_set = set(filenames)
    requested_list = sorted(requested_set)
    args = ['kpsewhich'] + requested_list
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()
    p_output = p.stdout.read()
    p_error = p.stderr.read()
    # get the subset of installed package names
    installed_set = set()
    for raw_line in p_output.splitlines():
        line = raw_line.strip()
        name = os.path.basename(line)
        if name in requested_set:
            installed_set.add(name)
        else:
            raise CheckPackageError(line)
    # return the set of installed filenames
    return installed_set

def check_installation(class_names, package_names):
    """
    @param class_names: a collection of class names
    @param package_names: a collection of package names
    @return: (installed_class_set, installed_package_set)
    """
    filenames = []
    filenames.extend(s + '.cls' for s in class_names)
    filenames.extend(s + '.sty' for s in package_names)
    installed_filenames = _check_installed_files(filenames)
    i_class = [s for s in class_names if s + '.cls' in installed_filenames]
    i_packs = [s for s in package_names if s + '.sty' in installed_filenames]
    return set(i_class), set(i_packs)

def check_packages(package_names):
    class_names = set([])
    installed_classes, installed_packages = check_installation(
            class_names, package_names)
    return installed_packages

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

def _create_temp_pdf_file(latex_text, output_format='pdf'):
    """
    The returned path name base does not yet have the pdf or dvi extension.
    @param latex_text: contents of a LaTeX file
    @param output_format: either pdf or dvi
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
            '-output-format=%s' % output_format,
            '-halt-on-error',
            pathname]
    # 
    proc = subprocess.Popen(args,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    proc_stdout, proc_stderr = proc.communicate()
    return pathname

def get_png_contents_imagemagick(latex_text):
    """
    This involves lots of temporary files.
    pdflatex foo.tex
    convert -density 300 foo.pdf foo.png
    """
    # make the pdf file
    pathname = _create_temp_pdf_file(latex_text)
    # make the png file
    cmd = [
        'convert',
        #'-density', '100',
        pathname + '.pdf',
        pathname + '.png']
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    proc_stdout, proc_stderr = proc.communicate()
    if proc.returncode:
        raise ValueError(proc_stderr)
    # read the png file
    with open(pathname + '.png', 'rb') as fin:
        image_data = fin.read()
    return image_data

def get_png_contents_ghostscript(latex_text):
    """
    This seemed to give ugly looking png files for some reason.
    It involves using temporary files.
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
    # sDEVICE used to be png16m
    args = [
            'gs', '-dUseCropBox', '-r144',
            '-dSAFER', '-dBATCH', '-dNOPAUSE', '-sDEVICE=pngalpha',
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
    except IOError as e:
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
    except IOError as e:
        raise ValueError('failed to create a pdf file')
    return pdf_contents

def latex_text_to_response(latex_text, latexformat):
    """
    @param latex_text: the text of a latex file
    @param latexformat: one of three possible formats
    @return: a response suitable to return from the get_response interface
    """
    assert_latexformat(latexformat)
    if latexformat == LATEXFORMAT_TEX:
        return latex_text
    elif latexformat == LATEXFORMAT_PDF:
        return get_pdf_contents(latex_text)
    elif latexformat == LATEXFORMAT_PNG:
        return get_png_contents_ghostscript(latex_text)

def get_response(
        requested_documentclass, document_body, latexformat,
        packages=(), preamble=''):
    """
    @param requested_documentclass: the documentclass
    @param document_body: the text inside a document environment
    @param latexformat: one of three latex output formats
    @param packages: a collection of requested packages
    @param preamble: color definitions, for example
    @return: a response suitable to return from the get_response interface
    """
    # check the requested format
    assert_latexformat(latexformat)
    # get the subset of installed class and package names
    requested_class_names = set([requested_documentclass])
    requested_package_names = set(packages)
    installed_class_names, installed_package_names = check_installation(
            requested_class_names, requested_package_names)
    # If no essential packages are missing
    # then define the documentclass and usepackage lines.
    missing_list = sorted(requested_package_names - installed_package_names)
    if not missing_list:
        if requested_documentclass in installed_class_names:
            documentclass = requested_documentclass
        else:
            documentclass = 'article'
        usepackage_list = sorted(installed_package_names)
    # If essential packages are missing and a compiled response
    # has been requested then this is a problem.
    # If essential packages are missing but a compiled response
    # has not been requested, then show the most optimistic settings.
    if missing_list:
        if latexformat in (LATEXFORMAT_PDF, LATEXFORMAT_PNG):
            msg = 'missing LaTeX packages: ' + ' '.join(missing_list)
            raise LatexPackageError(msg)
        else:
            documentclass = requested_documentclass
            usepackage_list = sorted(requested_set)
    # define the usepackage text
    usepackage_lines = ['\\usepackage{%s}' % s for s in usepackage_list]
    usepackage_text = '\n'.join(usepackage_lines)
    # define the latex body
    chunks = [
        '\\documentclass{%s}' % documentclass,
        usepackage_text,
        preamble,
        '\\begin{document}',
        document_body,
        '\\end{document}']
    latex_text = '\n'.join(c for c in chunks if c)
    # respond using the requested format
    return latex_text_to_response(latex_text, latexformat)

def get_centered_figure_response(
        figure_body, latexformat, figure_caption, figure_label,
        packages=(), preamble=''):
    """
    @param requested_documentclass: the documentclass
    @param figure_body: the contents of a centered figure environment
    @param latexformat: one of three latex output formats
    @param figure_caption: caption of the figure or None
    @param figure_label: label of the figure or None
    @param packages: a collection of requested packages
    @param preamble: color definitions, for example
    @return: a response suitable to return from the get_response interface
    """
    # TODO check the following assumption
    # When a figure is requested always use an article documentclass.
    documentclass = 'article'
    # The documentclass and usepackage list do not depend on the installation.
    # Therefore only check for missing packages if compilation is requested.
    requested_pkg_set = set(packages)
    installed_pkg_set = check_packages(requested_pkg_set)
    if latexformat in (LATEXFORMAT_PDF, LATEXFORMAT_PNG):
        missing_list = sorted(requested_pkg_set - installed_pkg_set)
        if missing_list:
            msg = 'missing LaTeX packages: ' + ' '.join(missing_list)
            raise LatexPackageError(msg)
    # define the usepackage text
    usepackage_list = sorted(requested_pkg_set)
    usepackage_lines = ['\\usepackage{%s}' % s for s in usepackage_list]
    usepackage_text = '\n'.join(usepackage_lines)
    # define the caption and label of the figure if they exist
    caption_chunk = ''
    if figure_caption:
        caption_chunk = '\\caption{\n' + figure_caption + '\n}'
    label_chunk = ''
    if figure_label:
        label_chunk = '\\label{' + figure_label + '}'
    # define the latex body
    chunks = [
        '\\documentclass{%s}' % documentclass,
        usepackage_text,
        preamble,
        '\\begin{document}',
        '\\begin{figure}',
        '\\centering',
        figure_body,
        caption_chunk,
        label_chunk,
        '\\end{figure}',
        '\\end{document}']
    latex_text = '\n'.join(c for c in chunks if c)
    # respond using the requested format
    return latex_text_to_response(latex_text, latexformat)


class TestLatexUtil(unittest.TestCase):

    def test_options_dict_to_string_nonempty(self):
        observed = options_to_string('auto', scale=0.5)
        expected = '[auto,scale=0.5]'
        self.assertEqual(observed, expected)

    def test_options_dict_to_string_empty(self):
        observed = options_to_string()
        expected = ''
        self.assertEqual(observed, expected)


if __name__ == '__main__':
    unittest.main()
