""" Draw a tree annotated with roots of eigenfunctions.

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
            Form.Integer('eig_idx1',
                'first eigenfunction index (1 means Fiedler)', 1, low=0),
            Form.Integer('eig_idx2',
                'second eigenfunction index (1 means Fiedler)', 2, low=0),
            Form.ImageFormat(),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.Image('tree')

def get_response_content(fs):
    # get a properly formatted newick tree with branch lengths
    tree = Newick.parse(fs.tree, SpatialTree.SpatialTree)
    # get the vertex valuations
    id_to_v1 = Harmonic.get_harmonic_valuations(tree, fs.eig_idx1)
    id_to_v2 = Harmonic.get_harmonic_valuations(tree, fs.eig_idx2)
    # do the layout
    try:
        layout = FastDaylightLayout.StraightBranchLayout()
        layout.do_layout(tree)
    except RuntimeError as e:
        pass
    # draw the image
    try:
        ext = Form.g_imageformat_to_ext[fs.imageformat]
        return DrawEigenLacing.get_single_tree_image(
                tree, (640, 480), ext, id_to_v1, id_to_v2)
    except CairoUtil.CairoUtilError as e:
        raise HandlingError(e)
