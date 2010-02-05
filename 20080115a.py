"""Given a newick tree, check it for syntax errors and briefly summarize it.
"""

from StringIO import StringIO

from SnippetUtil import HandlingError
import Newick
import Form

def get_form():
    """
    @return: the body of a form
    """
    # define the newick string
    tree_string = Newick.daylight_example_tree
    tree = Newick.parse(tree_string, Newick.NewickTree)
    formatted_tree_string = Newick.get_narrow_newick_string(tree, 60)
    # return the form objects
    return [Form.MultiLine('tree', 'newick tree', formatted_tree_string)]

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # get the tree
    tree = Newick.parse(fs.tree, Newick.NewickTree)
    tree.assert_valid()
    # write the response
    out = StringIO()
    print >> out, '\n'.join(tree.gen_description_lines())
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
