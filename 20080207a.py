"""Sample joint nucleotide substitutions on a fixed tree with known tips.

Maybe this is the one that was buggy.
"""

import numpy as np

from SnippetUtil import HandlingError
import SnippetUtil
import Newick
import SpatialTree
import DrawTreeImage
import EqualArcLayout
import PhyLikelihood
import RateMatrix
import MatrixUtil
import Util
import PathSampler
import Form
import FormOut

# This example tree might be useful for debugging.
# (((a:2, b:2)A:2)X:2, ((c:2)x:2, d:2)B:2);
# ((a:1, b:2)A:2, (c:3, d:4)B:2);
# ((a:1, b:2)A:1, (c:3, d:4)B:1, (e:0.25, f:0.5)C:1);
g_tree_string = '((a:1, b:2)A:1, (c:3, d:4)B:1, (e:0.25, f:0.5)C:1);'

def get_form():
    """
    @return: the body of a form
    """
    # define the tree string
    tree = Newick.parse(g_tree_string, Newick.NewickTree)
    formatted_tree_string = Newick.get_narrow_newick_string(tree, 60)
    # define the default tip data lines
    default_tip_data_lines = [
            'a : A',
            'b : A',
            'c : C',
            'd : T',
            'e : T',
            'f : T']
    # define the default rate matrix lines
    R = np.array([
        [-1, 1/3.0, 1/3.0, 1/3.0],
        [1/3.0, -1, 1/3.0, 1/3.0],
        [1/3.0, 1/3.0, -1, 1/3.0],
        [1/3.0, 1/3.0, 1/3.0, -1]])
    # define the form objects
    form_objects = [
            Form.MultiLine('tree', 'newick tree', formatted_tree_string),
            Form.MultiLine('column', 'tip data',
                '\n'.join(default_tip_data_lines)),
            Form.Matrix('matrix', 'rate matrix',
                R, MatrixUtil.assert_rate_matrix),
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
    # get the dictionary mapping the branch name to the nucleotide
    name_to_nt = {}
    lines = Util.get_stripped_lines(fs.column.splitlines())
    if lines:
        name_to_nt = SnippetUtil.get_generic_dictionary(lines, 'name',
                'nucleotide', list('acgtACGT'))
    # augment the tips with the nucleotide letters
    for name, nt in name_to_nt.items():
        try:
            node = tree.get_unique_node(name)
        except Newick.NewickSearchError as e:
            raise HandlingError(e)
        if node.children:
            msg = 'constraints on internal nodes are not implemented'
            raise HandlingError(msg)
        node.state = nt.upper()
    # read the rate matrix
    R = fs.matrix
    # convert the rate matrix to a rate matrix object
    states = list('ACGT')
    rate_matrix_object = RateMatrix.RateMatrix(R.tolist(), states)
    # simulate the ancestral nucleotides
    rate_matrix_object.simulate_ancestral_states(tree)
    # simulate a path on each branch
    # this breaks up the branch into a linear sequence of nodes and adds color
    for node in tree.gen_non_root_nodes():
        simulate_branch_path(tree, node, rate_matrix_object)
    # do the layout
    EqualArcLayout.do_layout(tree)
    # draw the image
    try:
        ext = Form.g_imageformat_to_ext[fs.imageformat]
        return DrawTreeImage.get_tree_image(tree, (640, 480), ext)
    except CairoUtil.CairoUtilError as e:
        raise HandlingError(e)

def simulate_branch_path(tree, node, rate_matrix_object):
    # purines are red; pyrimidines are blue
    # A and T are brighter, G and C are darker
    nt_to_color = {'A':'ff4444', 'G':'ffaaaa', 'T':'4444ff', 'C':'aaaaff'}
    node.branch_color = nt_to_color[node.state]
    rate_matrix = rate_matrix_object.dictionary_rate_matrix
    initial_state = node.parent.state
    terminal_state = node.state
    states = 'ACGT'
    events = None
    while events is None:
        events = PathSampler.get_nielsen_sample(
                initial_state, terminal_state, states, node.blen, rate_matrix)
    parent = node.parent
    last_t = 0
    for t, state in events:
        new = SpatialTree.SpatialTreeNode()
        new.name = node.name
        new.state = state
        new.branch_color = nt_to_color[parent.state]
        tree.insert_node(new, parent, node, (t - last_t) / float(node.blen))
        last_t = t
        parent = new
