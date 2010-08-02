"""Given a newick tree, draw an ascii representation.
"""

import Newick
import DrawTree
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
    form_objects = [
            Form.MultiLine('tree', 'newick tree', formatted_tree_string),
            Form.RadioGroup('group1', 'output options', [
                Form.RadioItem('choice1',
                    'a simple indented format', True),
                Form.RadioItem('choice2',
                    'ascii art with branch lengths'),
                Form.RadioItem('choice3',
                    'ascii art with topology only and less whitespace'),
                Form.RadioItem('choice4',
                    'ascii art with topology only and more whitespace')])]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):
    tree = Newick.parse(fs.tree, Newick.NewickTree)
    tree.assert_valid()
    # Draw the tree using a simple indented format
    # if the user chose the first option.
    # Otherwise create a drawing object and draw fancier ascii art.
    if fs.choice1:
        return str(tree) + '\n'
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
        return drawer.draw(tree) + '\n'
