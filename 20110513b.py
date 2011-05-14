"""
Draw trees annotated with roots of eigenfunctions, revised for TikZ.

Use the harmonic extensions of eigenvectors
of the Schur complement Laplacian matrix.
The pgf manual pgfmanualse13 has helpful information
about TikZ coordinate transformations.
"""

from SnippetUtil import HandlingError
import Newick
import SpatialTree
import FastDaylightLayout
import Form
import FormOut
import DrawEigenLacing
import Harmonic
import tikz


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
                Form.CheckItem('reflect_trees', 'reflect trees', True)]),
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

def get_response_content(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: the response
    """
    # get a properly formatted newick tree with branch lengths
    tree = Newick.parse(fs.tree, SpatialTree.SpatialTree)
    # check the indices
    if fs.last_index <= fs.first_index:
        msg = 'the last index should be greater than the first index'
        raise ValueError(msg)
    # get the vertex valuations
    valuations = [Harmonic.get_harmonic_valuations(
        tree, i) for i in range(fs.first_index, fs.last_index+1)]
    # do the layout
    try:
        layout = FastDaylightLayout.StraightBranchLayout()
        layout.do_layout(tree)
    except RuntimeError, e:
        pass
    # draw the image
    physical_size = (fs.width, fs.height)
    tikz_text = DrawEigenLacing.get_forest_image_tikz(
            tree, physical_size, valuations, fs.inner_margin, fs.reflect_trees)
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

