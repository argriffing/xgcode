"""Given a newick tree, draw an image using one of several file formats.
"""

from SnippetUtil import HandlingError
import Newick
import SpatialTree
import DrawTreeImage
import EqualArcLayout
import FastDaylightLayout
import CairoUtil
import Form
import FormOut

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
            Form.RadioGroup('layout', 'tree layout options', [
                Form.RadioItem('daylight', 'equal daylight layout', True),
                Form.RadioItem('arc', 'equal arc layout'),
                Form.RadioItem('curved', 'curved branch layout')]),
            Form.ImageFormat(),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.Image('tree')

def get_response_content(fs):
    # get a properly formatted newick tree with branch lengths
    tree = Newick.parse(fs.tree, SpatialTree.SpatialTree)
    tree.assert_valid()
    if tree.has_negative_branch_lengths():
        msg = 'drawing a tree with negative branch lengths is not implemented'
        raise HandlingError(msg)
    tree.add_branch_lengths()
    # do the layout
    if fs.daylight:
        try:
            layout = FastDaylightLayout.StraightBranchLayout()
            layout.do_layout(tree)
        except RuntimeError as e:
            pass
    elif fs.curved:
        try:
            layout = FastDaylightLayout.CurvedBranchLayout()
            layout.set_min_segment_count(400)
            layout.do_layout(tree)
        except RuntimeError as e:
            pass
    elif fs.arc:
        EqualArcLayout.do_layout(tree)
    # draw the image
    try:
        ext = Form.g_imageformat_to_ext[fs.imageformat]
        return DrawTreeImage.get_tree_image(tree, (640, 480), ext)
    except CairoUtil.CairoUtilError as e
        raise HandlingError(e)
