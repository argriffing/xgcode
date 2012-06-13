"""
Draw trees annotated with roots of eigenfunctions, revised.

Use the harmonic extensions of eigenvectors
of the Schur complement Laplacian matrix.
"""

from SnippetUtil import HandlingError
import Newick
import SpatialTree
import FastDaylightLayout
import Form
import FormOut
import CairoUtil
import DrawEigenLacing
import Harmonic


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
            Form.IntegerInterval(
                'first_index', 'last_index',
                'eigenfunction index range (1 means Fiedler)',
                0, 12, low=0),
            Form.CheckGroup('check_options', 'output options', [
                Form.CheckItem('reflect_trees', 'reflect trees', True),
                Form.CheckItem('draw_background', 'draw background', True),
                Form.CheckItem('draw_labels', 'draw labels', True)]),
            Form.Integer('width', 'physical width', 880, low=10, high=1000),
            Form.Integer('height', 'physical height', 680, low=10, high=1000),
            Form.Integer('inner_margin', 'inner margin', 10, low=0, high=1000),
            Form.Integer('outer_margin', 'outer margin', 0, low=0, high=1000),
            Form.ImageFormat()]
    return form_objects

def get_form_out():
    return FormOut.Image('tree')

def get_response_content(fs):
    # get a properly formatted newick tree with branch lengths
    tree = Newick.parse(fs.tree, SpatialTree.SpatialTree)
    # get the vertex valuations
    valuations = [Harmonic.get_harmonic_valuations(
        tree, i) for i in range(fs.first_index, fs.last_index+1)]
    # do the layout
    try:
        layout = FastDaylightLayout.StraightBranchLayout()
        layout.do_layout(tree)
    except RuntimeError as e:
        pass
    # draw the image
    try:
        ext = Form.g_imageformat_to_ext[fs.imageformat]
        physical_size = (fs.width, fs.height)
        return DrawEigenLacing.get_forest_image_revised(
                tree, physical_size, ext, valuations,
                fs.draw_background, fs.draw_labels,
                fs.inner_margin, fs.outer_margin, fs.reflect_trees)
    except CairoUtil.CairoUtilError as e:
        raise HandlingError(e)

