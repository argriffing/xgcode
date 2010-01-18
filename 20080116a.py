"""Given a newick tree, draw an ascii representation.
"""

import StringIO

from SnippetUtil import HandlingError
import Newick
import DrawTree
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
            Form.RadioGroup('group1', 'output options', [
                Form.RadioItem('choice1', 'a simple indented format', True),
                Form.RadioItem('choice2', 'ascii art with branch lengths'),
                Form.RadioItem('choice3', 'ascii art with topology only and less whitespace'),
                Form.RadioItem('choice4', 'ascii art with topology only and more whitespace')])]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # get the tree
    tree = Newick.parse(fs.tree, Newick.NewickTree)
    tree.assert_valid()
    # start writing the output
    out = StringIO.StringIO()
    # Draw the tree using a simple indented format if the user chose the first option.
    # Otherwise create a drawing object and draw fancier ascii art.
    if fs.choice1:
        print >> out, tree
    else:
        # customize the tree drawing object according to the user choice
        drawer = DrawTree.DrawTree()
        if fs.choice2:
            drawer.use_branch_lengths = True
            drawer.force_ultrametric = False
            drawer.vertical_spacing = 1
            drawer.horizontal_spacing = 1
        elif fs.choice3:
            drawer.use_branch_lengths = False
            drawer.force_ultrametric = False
            drawer.vertical_spacing = 1
            drawer.horizontal_spacing = 1
        elif fs.choice4:
            drawer.use_branch_lengths = False
            drawer.force_ultrametric = True
            drawer.vertical_spacing = 3
            drawer.horizontal_spacing = 6
        # draw the tree
        print >> out, drawer.draw(tree)
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().rstrip()
