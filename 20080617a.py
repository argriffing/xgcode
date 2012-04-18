"""Given a newick tree, draw an image highlighting selected taxa.
"""

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
import FormOut

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
            Form.ImageFormat()]
    return form_objects

def get_form_out():
    return FormOut.Image('tree')

def split_branches(tree):
    """
    Split the branches of a tree.
    @param tree: a newick tree
    """
    old_nodes = list(tree.preorder())
    for node in old_nodes:
        if node is tree.root:
            if node.blen is not None:
                msg = 'the root node should not have a branch length'
                raise HandlingError(msg)
        elif node.blen is None:
            msg = 'each non-root node should have a branch length'
            raise HandlingError(msg)
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

def get_response_content(fs):
    # start writing the response type
    response_headers = []
    # get a properly formatted newick tree with branch lengths
    tree = NewickIO.parse(fs.tree, SpatialTree.SpatialTree)
    tree.assert_valid()
    if tree.has_negative_branch_lengths():
        msg = 'drawing a tree with negative branch lengths is not implemented'
        raise HandlingError(msg)
    tree.add_branch_lengths()
    # get the selected taxa
    selected_taxa = Util.get_stripped_lines(fs.selection.splitlines())
    # verify that each name is present in the tree
    for name in selected_taxa:
        try:
            node = tree.get_unique_node(name)
        except NewickIO.NewickSearchError as e:
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
    except RuntimeError as e:
        pass
    # draw the image
    try:
        ext = Form.g_imageformat_to_ext[fs.imageformat]
        return DrawTreeImage.get_tree_image(tree, (640, 480), ext)
    except CairoUtil.CairoUtilError as e:
        raise HandlingError(e)
