""" Draw trees annotated with roots of eigenfunctions.

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
                'last eigenfunction index (1 means Fiedler)', 4, low=0),
            Form.Integer('ncols', 'use this many columns', 1, low=1),
            Form.CheckGroup('check_options', 'output options', [
                Form.CheckItem('draw_background', 'draw background', True)]),
            Form.ImageFormat(),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.Image('tree')

def get_response_content(fs):
    # get a properly formatted newick tree with branch lengths
    tree = Newick.parse(fs.tree, SpatialTree.SpatialTree)
    # check the indices
    if fs.last_index <= fs.first_index:
        msg = 'the last index should be greater than the first index'
        raise ValueError(msg)
    # get the vertex valuations
    valuations = [DrawEigenLacing.get_harmonic_valuations(
        tree, i) for i in range(fs.first_index, fs.last_index+1)]
    # do the layout
    try:
        layout = FastDaylightLayout.StraightBranchLayout()
        layout.do_layout(tree)
    except RuntimeError, e:
        pass
    # draw the image
    try:
        ext = Form.g_imageformat_to_ext[fs.imageformat]
        return DrawEigenLacing.get_forest_image(
                tree, (640, 480), ext, valuations,
                fs.ncols, fs.draw_background)
    except CairoUtil.CairoUtilError, e:
        raise HandlingError(e)
