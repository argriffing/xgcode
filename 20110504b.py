"""Make a TikZ tree that shows harmonically extended valuations.
"""

from StringIO import StringIO
import argparse

from SnippetUtil import HandlingError
import Form
import FormOut
import tikz
import Newick
import TreeProjection

g_default_yaw = 0
g_default_pitch = 0.2
g_default_eigenvector_index = 2

#g_default_tree = '((1:1, 2:0.5)6:1, (3:0.333333333333, 4:0.5)7:1, 5:1)8;' 
g_default_tree = Newick.daylight_example_tree

def get_form():
    """
    @return: the body of a form
    """
    # define the tree string
    tree = Newick.parse(g_default_tree, Newick.NewickTree)
    formatted_tree_string = Newick.get_narrow_newick_string(tree, 60)
    # get the form objects
    form_objects = [
            Form.MultiLine('tree', 'tree', formatted_tree_string),
            Form.Integer('eigenvector_index',
                'eigenvector index (1 is Fiedler)',
                g_default_eigenvector_index, low=1),
            Form.Float('yaw', 'yaw', g_default_yaw),
            Form.Float('pitch', 'pitch', g_default_pitch),
            #Form.Float('scaling_factor',
            #'tikz scaling factor', 1.0, low_exclusive=0),
            Form.TikzFormat(),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.Tikz()

def get_tikz_text(tikz_body, scaling_factor):
    if scaling_factor != 1:
        sf = ',scale=%s' % scaling_factor
    else:
        sf = ''
    #tikz_header = r'\begin{tikzpicture}[auto%s]' % sf
    tikz_header = r'\begin{tikzpicture}[x=1em,y=1em,inner sep=0pt%s]' % sf
    tikz_footer = r'\end{tikzpicture}'
    return '\n'.join([tikz_header, tikz_body, tikz_footer])

def get_latex_text(tikz_text):
    latex_header = '\n'.join([
        '\\documentclass{article}',
        '\\usepackage{tikz}',
        '\\begin{document}'])
    latex_body = tikz_text
    latex_footer = r'\end{document}'
    return '\n'.join([latex_header, latex_body, latex_footer])

def get_response_content(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: the response
    """
    # get the texts
    tikz_lines = TreeProjection.get_tikz_lines(
        fs.tree, fs.eigenvector_index, fs.yaw, fs.pitch)
    tikz_text = get_tikz_text('\n'.join(tikz_lines), 1.0)
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

def main(args):
    tikz_lines = TreeProjection.get_tikz_lines(
        args.tree, args.eigenvector_index,
        args.yaw, args.pitch)
    tikz_text = get_tikz_text('\n'.join(tikz_lines), 1.0)
    latex_text = get_latex_text(tikz_text)
    print latex_text

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tree', default=g_default_tree,
            help='newick tree with branch lengths')
    parser.add_argument('--eigenvector_index',
            default=g_default_eigenvector_index, type=int,
            help='eigenvector index (1 is Fiedler)')
    parser.add_argument('--yaw', default=g_default_yaw, type=float),
    parser.add_argument('--pitch', default=g_default_pitch, type=float),
    main(parser.parse_args())

