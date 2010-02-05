"""Given a newick tree, draw an image using one of several file formats.
"""

from StringIO import StringIO
import random

from SnippetUtil import HandlingError
import SnippetUtil
import Newick
import SpatialTree
import DrawTreeImage
import EqualArcLayout
import FastDaylightLayout
import CairoUtil
import Form

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
            Form.RadioGroup('imageformat', 'image format options', [
                Form.RadioItem('png', 'png', True),
                Form.RadioItem('svg', 'svg'),
                Form.RadioItem('pdf', 'pdf'),
                Form.RadioItem('ps', 'ps')]),
            Form.RadioGroup('contentdisposition', 'image delivery options', [
                Form.RadioItem('inline', 'view the image', True),
                Form.RadioItem('attachment', 'download the image')])]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # start writing the response type
    response_headers = []
    # get a properly formatted newick tree with branch lengths
    tree = Newick.parse(fs.tree, SpatialTree.SpatialTree)
    tree.assert_valid()
    if tree.has_negative_branch_lengths():
        raise HandlingError('drawing a tree with negative branch lengths is not implemented')
    tree.add_branch_lengths()
    # do the layout
    if fs.daylight:
        try:
            layout = FastDaylightLayout.StraightBranchLayout()
            layout.do_layout(tree)
        except RuntimeError, e:
            pass
    elif fs.curved:
        try:
            layout = FastDaylightLayout.CurvedBranchLayout()
            layout.set_min_segment_count(400)
            layout.do_layout(tree)
        except RuntimeError, e:
            pass
    elif fs.arc:
        EqualArcLayout.do_layout(tree)
    # draw the image
    try:
        image_string = DrawTreeImage.get_tree_image(tree, (640, 480), fs.imageformat)
    except CairoUtil.CairoUtilError, e:
        raise HandlingError(e)
    # specify the content type
    format_to_content_type = {'svg':'image/svg+xml', 'png':'image/png', 'pdf':'application/pdf', 'ps':'application/postscript'}
    response_headers.append(('Content-Type', format_to_content_type[fs.imageformat]))
    # specify the content disposition
    image_filename = 'tree.' + fs.imageformat
    response_headers.append(('Content-Disposition', "%s; filename=%s" % (fs.contentdisposition, image_filename)))
    # return the response
    return response_headers, image_string
