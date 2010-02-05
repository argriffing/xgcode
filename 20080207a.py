"""Sample joint nucleotide substitutions on a fixed tree with known tips.
"""

from StringIO import StringIO

import numpy

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
    R = numpy.array([
        [-1, 1/3.0, 1/3.0, 1/3.0],
        [1/3.0, -1, 1/3.0, 1/3.0],
        [1/3.0, 1/3.0, -1, 1/3.0],
        [1/3.0, 1/3.0, 1/3.0, -1]])
    # define the form objects
    form_objects = [
            Form.MultiLine('tree', 'newick tree', formatted_tree_string),
            Form.MultiLine('column', 'tip data', '\n'.join(default_tip_data_lines)),
            Form.Matrix('matrix', 'rate matrix', R, MatrixUtil.assert_rate_matrix),
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
    # get the dictionary mapping the branch name to the nucleotide
    name_to_nucleotide = {}
    lines = list(Util.stripped_lines(StringIO(fs.column)))
    if lines:
        name_to_nucleotide = SnippetUtil.get_generic_dictionary(lines, 'name', 'nucleotide', list('acgtACGT'))
    # augment the tips with the nucleotide letters
    for name, nucleotide in name_to_nucleotide.items():
        try:
            node = tree.get_unique_node(name)
        except Newick.NewickSearchError, e:
            raise HandlingError(e)
        if node.children:
            raise HandlingError('constraints on internal nodes are not implemented')
        node.state = nucleotide.upper()
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

def simulate_branch_path(tree, node, rate_matrix_object):
    try:
        import RateMatrix
        import PathSampler
        import SpatialTree
    except ImportError, e:
        raise HandlingError(e)
    # purines are red; pyrimidines are blue
    # A and T are brighter, G and C are darker
    nucleotide_to_color = {'A':'ff4444', 'G':'ffaaaa', 'T':'4444ff', 'C':'aaaaff'}
    node.branch_color = nucleotide_to_color[node.state]
    rate_matrix = rate_matrix_object.dictionary_rate_matrix
    initial_state = node.parent.state
    terminal_state = node.state
    states = 'ACGT'
    events = None
    while events is None:
        events = PathSampler.get_nielsen_sample(initial_state, terminal_state, states, node.blen, rate_matrix)
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
