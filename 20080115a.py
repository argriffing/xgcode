"""Given a newick tree, check it for syntax errors and briefly summarize it.
"""

import Newick
import Form
import FormOut

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

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):
    tree = Newick.parse(fs.tree, Newick.NewickTree)
    tree.assert_valid()
    return '\n'.join(tree.gen_description_lines()) + '\n'
