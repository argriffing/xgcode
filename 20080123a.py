""" Given a newick tree, draw an image using colored branches.
"""

from StringIO import StringIO

from SnippetUtil import HandlingError
import SnippetUtil
import Newick
import SpatialTree
import DrawTreeImage
import FastDaylightLayout
import Form
import FormOut
import iterutils

def get_form():
    """
    @return: the body of a form
    """
    # define the tree string
    tree_string = Newick.daylight_example_tree
    tree = Newick.parse(tree_string, Newick.NewickTree)
    formatted_tree_string = Newick.get_narrow_newick_string(tree, 60)
    # define the default color lines
    default_color_lines = [
            'A : ff0000',
            'B : 00ff00',
            'C : 0000ff']
    # define the form objects
    form_objects = [
            Form.MultiLine('tree', 'newick tree', formatted_tree_string),
            Form.MultiLine('coloration', 'branch colors',
                '\n'.join(default_color_lines)),
            Form.ImageFormat(),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.Image('tree', [])

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # get a properly formatted newick tree with branch lengths
    tree = Newick.parse(fs.tree, SpatialTree.SpatialTree)
    tree.assert_valid()
    if tree.has_negative_branch_lengths():
        msg = 'drawing a tree with negative branch lengths is not implemented'
        raise HandlingError(msg)
    tree.add_branch_lengths()
    # get the dictionary mapping the branch name to the rgb color
    name_to_rgb = {}
    # parse the coloration string
    for line in iterutils.stripped_lines(StringIO(fs.coloration)):
        # get the branch and its color
        name_string, rgb_string = SnippetUtil.get_state_value_pair(line)
        rgb_string = rgb_string.upper()
        # validate the rgb string
        if len(rgb_string) != 6:
            msg = 'expected each rgb string to be six characters long'
            raise HandlingError(msg)
        bad_letters = set(rgb_string) - set('0123456789ABCDEFabcdef')
        if bad_letters:
            msg = 'found invalid rgb characters: %s' % str(tuple(bad_letters))
            raise HandlingError(msg)
        # associate the branch with its color
        name_to_rgb[name_string] = rgb_string
    # color the branches
    for name, rgb in name_to_rgb.items():
        try:
            node = tree.get_unique_node(name)
        except Newick.NewickSearchError, e:
            raise HandlingError(e)
        node.branch_color = rgb
    # do the layout
    try:
        layout = FastDaylightLayout.StraightBranchLayout()
        layout.do_layout(tree)
    except RuntimeError, e:
        pass
    # get some options
    ext = Form.g_imageformat_to_ext[fs.imageformat]
    filename = 'tree.' + ext
    contenttype = Form.g_imageformat_to_contenttype[fs.imageformat]
    contentdisposition = '%s; filename=%s' % (fs.contentdisposition, filename)
    # draw the image
    try:
        image_string = DrawTreeImage.get_tree_image(tree, (640, 480), ext)
    except CairoUtil.CairoUtilError, e:
        raise HandlingError(e)
    response_headers = [
            ('Content-Type', contenttype),
            ('Content-Disposition', contentdisposition)]
    # return the response
    return response_headers, image_string
