"""Given a newick tree, break the tree into segments.
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
    tree_string = Newick.daylight_example_tree
    tree = Newick.parse(tree_string, Newick.NewickTree)
    formatted_tree_string = Newick.get_narrow_newick_string(tree, 60)
    # define the form objects
    form_objects = [
            Form.MultiLine('tree', 'newick tree', formatted_tree_string),
            Form.Integer('segments', 'minimum segment count',
                200, low=1, high=10000)]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # get the tree
    tree = Newick.parse(fs.tree, Newick.NewickTree)
    tree.assert_valid()
    # get the minimum number of segments
    min_segment_count = fs.segments
    # determine the maximum allowed branch length
    total_branch_length = tree.get_total_length()
    max_branch_length = total_branch_length / float(min_segment_count)
    # any branch longer than the max branch length will be broken in half
    while True:
        old_nodes = list(tree.preorder())
        for node in old_nodes:
            if node is tree.root:
                if node.blen is not None:
                    msg = 'the root node should not have a branch length'
                    raise HandlingError(msg)
            elif node.blen is None:
                msg = 'each non-root node should have a branch length'
                raise HandlingError(msg)
            elif node.blen > max_branch_length:
                # create a new node and set its attributes
                new = Newick.NewickNode()
                new.name = node.name
                # insert the new node
                tree.insert_node(new, node.parent, node, .5)
        # if no node was added then break out of the loop
        if len(old_nodes) == len(list(tree.preorder())):
            break
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, tree.get_newick_string()
