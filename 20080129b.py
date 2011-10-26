""" Given a tree and tip constraints, simulate JC nucleotide substitutions.

The tree should be in the Newick format.
JC is Jukes-Cantor.
This is almost certainly buggy and I think I knew the bug at one point.
The way it does the conditional probability is not right.
"""

# This example tree might be useful for debugging.
# (((a:2, b:2)A:2)X:2, ((c:2)x:2, d:2)B:2);
# ((a:1, b:2)A:2, (c:3, d:4)B:2);
# ((a:1, b:2)A:1, (c:3, d:4)B:1, (e:0.25, f:0.5)C:1);
g_tree_string = '((a:1, b:2)A:1, (c:3, d:4)B:1, (e:0.25, f:0.5)C:1);'

from SnippetUtil import HandlingError
import SnippetUtil
import RateMatrix
import PathSampler
import Newick
import SpatialTree
import DrawTreeImage
import EqualArcLayout
import PhyLikelihood
import MatrixUtil
import Form
import FormOut
import iterutils

def get_form():
    """
    @return: the body of a form
    """
    # define the tree string
    tree_string = g_tree_string
    tree = Newick.parse(tree_string, Newick.NewickTree)
    formatted_tree_string = Newick.get_narrow_newick_string(tree, 60)
    # define the default tip data lines
    default_tip_data_lines = [
            'a : A',
            'b : A',
            'c : C',
            'd : T',
            'e : T',
            'f : T']
    # define the form objects
    form_objects = [
            Form.MultiLine('tree', 'newick tree',
                formatted_tree_string),
            Form.MultiLine('column', 'tip data',
                '\n'.join(default_tip_data_lines)),
            Form.ImageFormat(),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.Image('tree')

def get_response(fs):
    # get a properly formatted newick tree with branch lengths
    tree = Newick.parse(fs.tree, SpatialTree.SpatialTree)
    tree.assert_valid()
    if tree.has_negative_branch_lengths():
        msg = 'drawing a tree with negative branch lengths is not implemented'
        raise HandlingError(msg)
    tree.add_branch_lengths()
    # get the dictionary mapping the branch name to the nucleotide
    name_to_nucleotide = {}
    # parse the column string
    for line in iterutils.stripped_lines(fs.column.splitlines()):
        name_string, nucleotide_string = SnippetUtil.get_state_value_pair(line)
        if nucleotide_string not in list('acgtACGT'):
            msg = '"%s" is not a valid nucleotide' % nucleotide_string
            raise HandlingError(msg)
        nucleotide_string = nucleotide_string.upper()
        if name_string in name_to_nucleotide:
            raise HandlingError('the name "%s" was duplicated' % name_string)
        name_to_nucleotide[name_string] = nucleotide_string
    # augment the tips with the nucleotide letters
    for name, nucleotide in name_to_nucleotide.items():
        try:
            node = tree.get_unique_node(name)
        except Newick.NewickSearchError as e:
            raise HandlingError(e)
        if node.children:
            msg = 'constraints on internal nodes are not implemented'
            raise HandlingError(msg)
        node.state = nucleotide
    # get the Jukes-Cantor rate matrix object
    dictionary_rate_matrix = RateMatrix.get_jukes_cantor_rate_matrix()
    ordered_states = list('ACGT')
    row_major_rate_matrix = MatrixUtil.dict_to_row_major(
            dictionary_rate_matrix, ordered_states, ordered_states)
    rate_matrix_object = RateMatrix.RateMatrix(
            row_major_rate_matrix, ordered_states)
    # simulate the ancestral nucleotides
    rate_matrix_object.simulate_ancestral_states(tree)
    # simulate a path on each branch
    # this breaks up the branch into a linear sequence of nodes and adds color
    for node in tree.gen_non_root_nodes():
        simulate_branch_path(tree, node)
    # do the layout
    EqualArcLayout.do_layout(tree)
    # draw the image
    try:
        ext = Form.g_imageformat_to_ext[fs.imageformat]
        return DrawTreeImage.get_tree_image(tree, (640, 480), ext)
    except CairoUtil.CairoUtilError as e:
        raise HandlingError(e)

def simulate_branch_path(tree, node):
    """
    Simulate the nucleotide history on the path between a node and its parent.
    This simulated path is conditional on known values at each node.
    Purines are red; pyrimidines are blue.
    A and T are brighter; G and C are darker.
    @param tree: a SpatialTree with simulated nucleotides at each node
    @param node: the node that defines the branch on which to simulate a history
    """
    nucleotide_to_color = {
            'A':'FF4444', 'G':'FF8888', 'T':'4444FF', 'C':'8888FF'}
    node.branch_color = nucleotide_to_color[node.state]
    rate_matrix = RateMatrix.get_jukes_cantor_rate_matrix()
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
        new.branch_color = nucleotide_to_color[parent.state]
        tree.insert_node(new, parent, node, (t - last_t) / float(node.blen))
        last_t = t
        parent = new
