"""Given a newick tree, calculate tip weights.
"""

from SnippetUtil import HandlingError
import Newick
import LeafWeights
import Form
import FormOut

def get_form():
    """
    @return: the body of a form
    """
    # define the tree string
    tree_string = LeafWeights.stone_example_tree
    tree = Newick.parse(tree_string, Newick.NewickTree)
    formatted_tree_string = Newick.get_narrow_newick_string(tree, 60)
    # define the form objects
    form_objects = [
            Form.MultiLine('tree', 'newick tree', formatted_tree_string),
            Form.RadioGroup('method', 'weight calculation method', [
                Form.RadioItem('stone',
                    'use the method of Stone and Sidow', True),
                Form.RadioItem('thompson',
                    'use the Thompson et al.')])]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):
    tree = Newick.parse(fs.tree, Newick.NewickTree)
    tree.assert_valid()
    tree.add_branch_lengths()
    if tree.has_negative_branch_lengths():
        msg_a = 'calculating weights for a tree '
        msg_b = 'with negative branch lengths is not implemented'
        raise HandlingError(msg_a + msg_b)
    if fs.stone:
        name_weight_pairs = LeafWeights.get_stone_weights(tree)
    elif fs.thompson:
        name_weight_pairs = LeafWeights.get_thompson_weights(tree)
    lines = ['%s: %f' % pair for pair in name_weight_pairs]
    return '\n'.join(lines) + '\n'
