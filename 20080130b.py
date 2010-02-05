"""Given a newick tree, merge segmented branches.
"""

from StringIO import StringIO

from SnippetUtil import HandlingError
import Newick
import Form

def get_form():
    """
    @return: the body of a form
    """
    # define the tree string
    tree_string = '((((a:5)a:5, (b:5)b:5)A:1.5)A:1.5, (((c:5)c:5, (d:5)d:5)B:1.5)B:1.5);'
    tree = Newick.parse(tree_string, Newick.NewickTree)
    formatted_tree_string = Newick.get_narrow_newick_string(tree, 60)
    # define the form objects
    return [Form.MultiLine('tree', 'newick tree', formatted_tree_string)]

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # get the tree
    tree = Newick.parse(fs.tree, Newick.NewickTree)
    tree.assert_valid()
    # modify the tree
    segmenting_nodes = [node for node in tree.preorder() if len(node.children) == 1]
    for node in segmenting_nodes:
        tree.remove_node(node)
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, tree.get_newick_string()
