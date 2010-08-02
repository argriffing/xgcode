"""Given a newick tree, remove branches with zero length.

Taxa with zero branch length to the rest of the tree may be removed.
"""

from SnippetUtil import HandlingError
import Newick
import Form
import FormOut

def get_form():
    """
    @return: the body of a form
    """
    default_tree_string = '((a:1, b:1):0, c:1);'
    return [Form.MultiLine('tree', 'newick tree', default_tree_string)]

def get_form_out():
    return FormOut.Newick()

def get_response_content(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # get the tree
    tree = Newick.parse(fs.tree, Newick.NewickTree)
    tree.assert_valid()
    # find nodes with a branch length of zero
    zero_blen_nodes = [node for node in tree.preorder() if node.blen == 0.0]
    # deal with these nodes appropriately
    for node in zero_blen_nodes:
        if node.blen != 0.0:
            # it is possible that this node no longer has branch length zero
            continue
        elif node.has_children():
            # remove internal nodes with zero branch length
            tree.remove_node(node)
        else:
            # prune terminal nodes with zero branch length
            intersection = tree.prune(node)
            # remove intersection nodes that have a single child
            if intersection.get_child_count() == 1:
                tree.remove_node(intersection)
    # return the response
    return tree.get_newick_string() + '\n'
