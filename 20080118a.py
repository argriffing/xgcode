"""Given a newick tree, calculate tip weights.
"""

from StringIO import StringIO

from SnippetUtil import HandlingError
import Newick
import LeafWeights
import Form

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
                Form.RadioItem('stone', 'use the method of Stone and Sidow', True),
                Form.RadioItem('thompson', 'use the Thompson et al.')])]
    return form_objects

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
        raise HandlingError('calculating weights for a tree with negative branch lengths is not implemented')
    # get the weights
    if fs.stone:
        name_weight_pairs = LeafWeights.get_stone_weights(tree)
    elif fs.thompson:
        name_weight_pairs = LeafWeights.get_thompson_weights(tree)
    # report the weights
    out = StringIO()
    for name, weight in name_weight_pairs:
        print >> out, '%s: %f' % (name, weight)
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
