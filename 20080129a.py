""" Given a newick tree, split each branch into two branches.
"""

from StringIO import StringIO

from SnippetUtil import HandlingError
import Newick
import Form
import FormOut

def get_form():
    """
    @return: the body of a form
    """
    # define the tree string
    tree_string = Newick.daylight_example_tree
    tree = Newick.parse(tree_string, Newick.NewickTree)
    formatted_tree_string = Newick.get_narrow_newick_string(tree, 60)
    # define the form objects
    return [Form.MultiLine('tree', 'newick tree', formatted_tree_string)]

def get_form_out():
    return FormOut.Newick()

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # get the tree
    tree = Newick.parse(fs.tree, Newick.NewickTree)
    tree.assert_valid()
    # modify the tree
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
            new = Newick.NewickNode()
            new.name = node.name
            # insert the new node
            tree.insert_node(new, node.parent, node, .5)
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, tree.get_newick_string()
