"""Given a newick tree, remove a node.
"""

from SnippetUtil import HandlingError
from LeafWeights import stone_example_tree
import Newick
import Form
import FormOut

def get_form():
    """
    @return: the body of a form
    """
    # define the tree string
    tree_string = stone_example_tree
    tree = Newick.parse(tree_string, Newick.NewickTree)
    formatted_tree_string = Newick.get_narrow_newick_string(tree, 60)
    # define the form objects
    form_objects = [
            Form.MultiLine('tree', 'newick tree', formatted_tree_string),
            Form.SingleLine('node', 'the name of the node to remove', 'a')]
    return form_objects

def get_form_out():
    return FormOut.Newick()

def get_response_content(fs):
    # get the tree
    tree = Newick.parse(fs.tree, Newick.NewickTree)
    tree.assert_valid()
    # get the node
    try:
        node = tree.get_unique_node(fs.node)
    except Newick.NewickSearchError, e:
        raise HandlingError(e)
    if node is tree.root:
        raise HandlingError('the root cannot be removed')
    # remove the node
    tree.remove_node(node)
    # return the response
    return tree.get_newick_string() + '\n'
