"""
Draw trees annotated with roots of eigenfunctions, revised for TikZ.

Use the harmonic extensions of eigenvectors
of the Schur complement Laplacian matrix.
The pgf manual pgfmanualse13 has helpful information
about TikZ coordinate transformations.
"""

import scipy
from scipy import linalg
import numpy as np

from SnippetUtil import HandlingError
import Newick
import SpatialTree
import FastDaylightLayout
import Form
import FormOut
import DrawEigenLacing
import tikz
import Ftree
import FtreeIO
import FtreeAux


def get_form():
    """
    @return: the body of a form
    """
    # define the tree string
    tree_string = Newick.daylight_example_tree
    tree = Newick.parse(tree_string, Newick.NewickTree)
    formatted_tree_string = Newick.get_narrow_newick_string(tree, 60)
    # define the form objects
    form_objects = [
            Form.MultiLine('tree', 'newick tree', formatted_tree_string),
            Form.Integer('first_index',
                'first eigenfunction index (1 means Fiedler)', 0, low=0),
            Form.Integer('last_index',
                'last eigenfunction index (1 means Fiedler)', 12, low=0),
            Form.CheckGroup('check_options', 'output options', [
                Form.CheckItem('reflect_trees',
                    'reflect trees across the vertical axis'),
                Form.CheckItem('show_subfigure_labels',
                    'show subfigure labels', True),
                Form.CheckItem('show_vertex_labels',
                    'show vertex labels')]),
            Form.Float('width', 'width (centimeters)',
                12, low_inclusive=1, high_exclusive=1000),
            Form.Float('height', 'height (centimeters)',
                8, low_inclusive=1, high_exclusive=1000),
            Form.Float('inner_margin', 'inner margin (centimeters)',
                0.25, low_inclusive=0, high_exclusive=1000),
            Form.TikzFormat(),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.Tikz('tkz')

def get_latex_text(latex_body):
    return '\n'.join([
        '\\documentclass{article}',
        '\\usepackage{tikz}',
        '\\usetikzlibrary{snakes}',
        '\\begin{document}',
        latex_body,
        '\\end{document}'])

def get_figure_text(figure_body):
    return '\n'.join([
        '\\begin{figure}',
        '\\centering',
        figure_body,
        '\\caption{',
        'This figure shows the interlacing',
        'of harmonically extended eigenvectors.',
        '}',
        '\\label{fig:interlacing}',
        '\\end{figure}'])

#FIXME this is used because the analogous Ftree function does not
#      include the constant vector
def TB_to_harmonic_valuations(T, B):
    """
    @param T: topology
    @param B: branch lengths
    @return: a number of dictionaries equal to the number of leaves
    """
    leaves = Ftree.T_to_leaves(T)
    internal = Ftree.T_to_internal_vertices(T)
    vertices = leaves + internal
    nleaves = len(leaves)
    Lbb = Ftree.TB_to_L_block(T, B, internal, internal)
    Lba = Ftree.TB_to_L_block(T, B, internal, leaves)
    L_schur = Ftree.TB_to_L_schur(T, B, leaves)
    w, v1 = scipy.linalg.eigh(L_schur)
    v2 = -np.dot(np.dot(np.linalg.pinv(Lbb), Lba), v1)
    V = np.vstack([v1, v2])
    vs = []
    for col in range(nleaves):
        d = dict((v, V[row, col]) for row, v in enumerate(vertices))
        vs.append(d)
    return vs

def get_response_content(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: the response
    """
    # get a properly formatted newick tree with branch lengths
    T, B, N = FtreeIO.newick_to_TBN(fs.tree)
    # check the indices
    if fs.last_index < fs.first_index:
        msg = 'the last index should not be greater than the first index'
        raise ValueError(msg)
    # get the vertex valuations
    all_valuations = TB_to_harmonic_valuations(T, B)
    valuations = all_valuations[fs.first_index:]
    nfigures = (fs.last_index - fs.first_index) + 1
    # do the layout
    v_to_location = FtreeAux.equal_daylight_layout(T, B, 3)
    # draw the image
    physical_size = (fs.width, fs.height)
    tikz_text = DrawEigenLacing.get_forest_image_ftree(
            T, B, N, v_to_location,
            physical_size, valuations, nfigures, fs.inner_margin,
            fs.reflect_trees, fs.show_vertex_labels, fs.show_subfigure_labels)
    latex_text = get_latex_text(get_figure_text(tikz_text))
    # decide the output format
    if fs.tikz:
        return tikz_text
    elif fs.tex:
        return latex_text
    elif fs.pdf:
        return tikz.get_pdf_contents(latex_text)
    elif fs.png:
        return tikz.get_png_contents(latex_text)

