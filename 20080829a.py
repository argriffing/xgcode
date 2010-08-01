"""Given a newick tree, calculate weights for a subset of the tips.
"""

from SnippetUtil import HandlingError
import Newick
import LeafWeights
import Util
from Form import RadioItem
import Form
import FormOut
import const

g_ucsc_tree_string = const.read('20100730x')

def get_form():
    """
    @return: the body of a form
    """
    # define the default tree string
    tree_string = g_ucsc_tree_string
    tree = Newick.parse(tree_string, Newick.NewickTree)
    formatted_tree_string = Newick.get_narrow_newick_string(tree, 60)
    # default name selection
    default_name_selection = ('human_hg18', 'chimp_panTro1', 'rabbit_oryCun1')
    # define the form objects
    form_objects = [
            Form.MultiLine('tree', 'newick tree', formatted_tree_string),
            Form.MultiLine('selection', 'selected tip names',
                '\n'.join(default_name_selection)),
            Form.RadioGroup('method', 'weighting method', [
                RadioItem('stone', 'use the method of Stone and Sidow', True),
                RadioItem('thompson', 'use the method of Thompson et al.')])]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # get the tree
    tree = Newick.parse(fs.tree, Newick.NewickTree)
    tree.assert_valid()
    tree.add_branch_lengths()
    if tree.has_negative_branch_lengths():
        msg_a = 'calculating weights for a tree '
        msg_b = 'with negative branch lengths is not implemented'
        raise HandlingError(msg_a + msg_b)
    # get the selected names
    selection = Util.get_stripped_lines(fs.selection.splitlines())
    selected_name_set = set(selection)
    possible_name_set = set(node.get_name() for node in tree.gen_tips())
    extra_names = selected_name_set - possible_name_set
    if extra_names:
        msg_a = 'the following selected names are not valid tips: '
        msg_b = str(tuple(extra_names))
        raise HandlingError(msg_a + msg_b)
    # prune the tree 
    for name in set(node.name for node in tree.gen_tips()) - set(selection): 
        try: 
            node = tree.get_unique_node(name) 
        except NewickSearchError, e: 
            raise HandlingError(e) 
        tree.prune(node)
    # get the weights
    if fs.stone:
        name_weight_pairs = LeafWeights.get_stone_weights(tree)
    elif fs.thompson:
        name_weight_pairs = LeafWeights.get_thompson_weights(tree)
    # report the weights
    lines = ['%s: %f' % pair for pair in name_weight_pairs]
    text = '\n'.join(lines) + '\n'
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, text
