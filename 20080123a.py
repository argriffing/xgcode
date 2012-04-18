""" Given a newick tree, draw an image using colored branches.
"""

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
            Form.ImageFormat()]
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
    # get the dictionary mapping the branch name to the rgb color
    name_to_rgb = {}
    # parse the coloration string
    for line in iterutils.stripped_lines(fs.coloration.splitlines()):
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
        except Newick.NewickSearchError as e:
            raise HandlingError(e)
        node.branch_color = rgb
    # do the layout
    try:
        layout = FastDaylightLayout.StraightBranchLayout()
        layout.do_layout(tree)
    except RuntimeError as e:
        pass
    # draw the image
    try:
        ext = Form.g_imageformat_to_ext[fs.imageformat]
        return DrawTreeImage.get_tree_image(tree, (640, 480), ext)
    except CairoUtil.CairoUtilError as e:
        raise HandlingError(e)
