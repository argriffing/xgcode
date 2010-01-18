"""Given a newick tree, scale the branch lengths.
"""

import cgi
import StringIO

from SnippetUtil import HandlingError
import SnippetUtil
import Newick
import Form

def get_form():
    """
    @return: the body of a form
    """
    # define the tree string
    tree_string = '(((Human:0.1, Chimpanzee:0.2):0.8, Gorilla:0.3):0.7, Orangutan:0.4, Gibbon:0.5);'
    tree = Newick.parse(tree_string, Newick.NewickTree)
    formatted_tree_string = Newick.get_narrow_newick_string(tree, 60)
    # define the form objects
    form_objects = [
            Form.MultiLine('tree', 'newick tree', formatted_tree_string),
            Form.Float('sf', 'scaling factor', 0.333333, low_exclusive=0)]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # get the newick string
    try:
        tree = Newick.parse(fs.tree, Newick.NewickTree)
        tree.assert_valid()
    except Newick.NewickSyntaxError, e:
        raise HandlingError(e)
    # scale the branch lengths
    for node in tree.preorder():
        if node.blen is not None:
            node.blen *= fs.sf
    # write the newick tree
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, tree.get_newick_string()
