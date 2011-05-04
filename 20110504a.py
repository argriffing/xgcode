"""Make more TikZ tree figures.
"""

from StringIO import StringIO

from SnippetUtil import HandlingError
from Form import RadioItem
import Form
import FormOut
import tikz

g_demo_lines = [
        r'\node (n217598284)[draw,shape=circle] at (-5, -3) {1};',
        r'\node (n217598252)[draw,shape=circle] at (-5, -1) {2};',
        r'\node (n217598380)[draw,shape=circle] at (-5, 1) {3};',
        r'\node (n217598828)[draw,shape=circle] at (-5, 3) {4};',
        r'\node (n217598188)[draw,shape=circle] at (-3, -2) {8};',
        r'\node (n217597740)[draw,shape=circle] at (-3, 2) {9};',
        r'\node (n217598348)[draw,shape=circle] at (-1, 0) {10};',
        r'\node (n217599116)[draw,shape=circle] at (6, 3) {5};',
        r'\node (n217599148)[draw,shape=circle] at (6, 1) {6};',
        r'\node (n217599180)[draw,shape=circle] at (4, -1) {7};',
        r'\node (n217599052)[draw,shape=circle] at (4, 2) {11};',
        r'\node (n217598988)[draw,shape=circle] at (2, 0) {12};',
        r'\path (n217597740) edge node {2} (n217598348);',
        r'\path (n217598828) edge node {2} (n217597740);',
        r'\path (n217597740) edge node {1} (n217598380);',
        r'\path (n217598252) edge node {2} (n217598188);',
        r'\path (n217598348) edge node {3} (n217598188);',
        r'\path (n217598188) edge node {8} (n217598284);',
        r'\path (n217598988) edge node {5} (n217599052);',
        r'\path (n217599052) edge node {5} (n217599116);',
        r'\path (n217599148) edge node {3} (n217599052);',
        r'\path (n217599180) edge node {2} (n217598988);',
        r'\path (n217598348) edge node {1} (n217598988);']


def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.Float('scaling_factor', 'scaling factor',
                1.0, low_exclusive=0),
            Form.TikzFormat(),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.Tikz()

def get_tikz_text(scaling_factor):
    if scaling_factor != 1:
        sf = ',scale=%s' % scaling_factor
    else:
        sf = ''
    tikz_header = r'\begin{tikzpicture}[auto%s]' % sf
    tikz_footer = r'\end{tikzpicture}'
    tikz_body = '\n'.join(g_demo_lines)
    return '\n'.join([tikz_header, tikz_body, tikz_footer])

def get_latex_text(scaling_factor):
    latex_header = '\n'.join([
        r'\documentclass{article}',
        r'\usepackage{tikz}',
        r'\begin{document}'])
    latex_footer = r'\end{document}'
    latex_body = get_tikz_text(scaling_factor)
    return '\n'.join([latex_header, latex_body, latex_footer])

def get_response_content(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: the response
    """
    # get the texts
    tikz_text = get_tikz_text(fs.scaling_factor)
    latex_text = get_latex_text(fs.scaling_factor)
    # decide the output format
    if fs.tikz:
        return tikz_text
    elif fs.tex:
        return latex_text
    elif fs.pdf:
        return tikz.get_pdf_contents(latex_text)
    elif fs.png:
        return tikz.get_png_contents(latex_text)

def main(options):
    print get_latex_text(0.5)

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    options, args = parser.parse_args()
    main(options)

