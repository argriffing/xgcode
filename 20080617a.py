"""Given a newick tree, draw an image highlighting selected taxa.
"""

import StringIO

from SnippetUtil import HandlingError
import Util
import SpatialTree
import DrawTreeImage
import FastDaylightLayout
import RateMatrix
import HeatMap
import FelTree
import NewickIO
import Form

def get_form():
    """
    @return: the body of a form
    """
    # define the default tree string
    tree_string = NewickIO.daylight_example_tree
    tree = NewickIO.parse(tree_string, FelTree.NewickTree)
    formatted_tree_string = NewickIO.get_narrow_newick_string(tree, 60)
    # define the form objects
    form_objects = [
            Form.MultiLine('tree', 'newick tree', formatted_tree_string),
            Form.MultiLine('selection', 'selected taxa', '\n'.join('ABFG')),
            Form.RadioGroup('imageformat', 'image format options', [
                Form.RadioItem('png', 'png', True),
                Form.RadioItem('svg', 'svg'),
                Form.RadioItem('pdf', 'pdf'),
                Form.RadioItem('ps', 'ps')]),
            Form.RadioGroup('contentdisposition', 'image delivery options', [
                Form.RadioItem('inline', 'view the image', True),
                Form.RadioItem('attachment', 'download the image')])]
    return form_objects

def split_branches(tree):
    """
    Split the branches of a tree.
    @param tree: a newick tree
    """
    old_nodes = list(tree.preorder())
    for node in old_nodes:
        if node is tree.root:
            if node.blen is not None:
                raise HandlingError('the root node should not have a branch length')
        elif node.blen is None:
            raise HandlingError('each non-root node should have a branch length')
        else:
            # create a new node and set its attributes
            new = SpatialTree.SpatialTreeNode()
            new.name = node.name
            # insert the new node
            tree.insert_node(new, node.parent, node, .5)

def add_colors(tree, selection):
    """
    Add branch colors to a newick tree.
    @param tree: a newick tree
    @param selection: a list of taxon names
    """
    # set the tip states
    for node in tree.gen_tips():
        if node.name in selection:
            node.state = 'a'
        else:
            node.state = 'b'
    # get the total length of the tree
    total_length = sum(node.blen for node in tree.gen_non_root_nodes())
    # define the rate matrix
    states = ('a', 'b')
    mu = 1.0 / total_length
    row_major_rate_matrix = [[-mu, mu], [mu, -mu]]
    rate_matrix_object = RateMatrix.RateMatrix(row_major_rate_matrix, states)
    # repeatedly reroot and calculate root state distributions
    internal_nodes = list(tree.gen_internal_nodes())
    for node in internal_nodes:
        tree.reroot(node)
        rate_matrix_object.add_probabilities(tree)
        weights = [node.state_to_subtree_prob[state] for state in states]
        node.state_distribution = Util.weights_to_distribution(weights)
    for node in tree.gen_tips():
        node.state_distribution = []
        for state in states:
            if state == node.state:
                node.state_distribution.append(1.0)
            else:
                node.state_distribution.append(0.0)
    # set the color of each branch
    for node in tree.gen_non_root_nodes():
        parent_probability = node.parent.state_distribution[0]
        current_probability = node.state_distribution[0]
        p = (parent_probability + current_probability) / 2.0
        r, g, b = HeatMap.blue_red_gradient(p)
        rgb_string = ('%02x%02x%02x' % (r, g, b)).upper()
        node.branch_color = rgb_string

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # start writing the response type
    response_headers = []
    # get a properly formatted newick tree with branch lengths
    tree = NewickIO.parse(fs.tree, SpatialTree.SpatialTree)
    tree.assert_valid()
    if tree.has_negative_branch_lengths():
        raise HandlingError('drawing a tree with negative branch lengths is not implemented')
    tree.add_branch_lengths()
    # get the selected taxa
    selected_taxa = list(Util.stripped_lines(StringIO.StringIO(fs.selection)))
    # verify that each name is present in the tree
    for name in selected_taxa:
        try:
            node = tree.get_unique_node(name)
        except NewickIO.NewickSearchError, e:
            raise HandlingError(e)
        if node.children:
            raise HandlingError('an internal node was selected: ' + name)
    # split the branches of the tree a couple of times
    split_branches(tree)
    split_branches(tree)
    # color the tree
    add_colors(tree, selected_taxa)
    # do the layout
    try:
        layout = FastDaylightLayout.StraightBranchLayout()
        layout.do_layout(tree)
    except RuntimeError, e:
        pass
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
