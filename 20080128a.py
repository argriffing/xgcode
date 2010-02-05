"""Given a newick tree, remove a set of tips.
"""

from StringIO import StringIO

from SnippetUtil import HandlingError
import Newick
import Util
import Form

def get_form():
    """
    @return: the body of a form
    """
    # define the tree string
    tree_string = '((a:10, b:10)A:3, (c:10, d:10)B:3);'
    # define the form objects
    form_objects = [
            Form.MultiLine('tree', 'newick tree', tree_string),
            Form.MultiLine('names', 'selected taxa', '\n'.join(('c', 'd'))),
            Form.RadioGroup('selection_style', 'selection style', [
                Form.RadioItem('remove', 'remove the selected taxa', True),
                Form.RadioItem('keep', 'keep the selected taxa')])]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # get the tree
    tree = Newick.parse(fs.tree, Newick.NewickTree)
    tree.assert_valid()
    # get the set of names
    selection = Util.get_stripped_lines(StringIO(fs.names))
    selected_name_set = set(selection)
    possible_name_set = set(node.get_name() for node in tree.gen_tips())
    extra_names = selected_name_set - possible_name_set
    if extra_names:
        raise HandlingError('the following selected names are not valid tips: %s' % str(tuple(extra_names)))
    # get the list of tip nodes to remove
    if fs.remove:
        nodes_to_remove = [node for node in tree.gen_tips() if node.name in selection]
    elif fs.keep:
        nodes_to_remove = [node for node in tree.gen_tips() if node.name not in selection]
    # prune the tree
    for node in nodes_to_remove:
        tree.prune(node)
    # merge segmented branches
    internal_nodes_to_remove = [node for node in tree.preorder() if node.get_child_count() == 1]
    for node in internal_nodes_to_remove:
        tree.remove_node(node)
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, tree.get_newick_string()
