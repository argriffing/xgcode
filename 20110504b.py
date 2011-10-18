"""Make a TikZ tree that shows harmonically extended valuations.
"""

from StringIO import StringIO

import Form
import FormOut
import tikz
import Newick
import TreeProjection

g_default_yaw = 0
g_default_pitch = 0.2
g_default_eigenvector_index = 2

g_default_tree = Newick.daylight_example_tree

def get_form():
    """
    @return: the body of a form
    """
    # define the tree string
    tree = Newick.parse(g_default_tree, Newick.NewickTree)
    formatted_tree_string = Newick.get_narrow_newick_string(tree, 60)
    # get the form objects
    form_objects = [
            Form.MultiLine('tree', 'tree', formatted_tree_string),
            Form.Integer('eigenvector_index',
                'eigenvector index (1 is Fiedler)',
                g_default_eigenvector_index, low=1),
            Form.Float('yaw', 'yaw', g_default_yaw),
            Form.Float('pitch', 'pitch', g_default_pitch),
            Form.TikzFormat(),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.Tikz()

def get_response_content(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: the response
    """
    tikz_body_lines = TreeProjection.get_tikz_lines(
            fs.tree, fs.eigenvector_index, fs.yaw, fs.pitch)
    tikz_body = '\n'.join(tikz_body_lines)
    options = {
            'x' : '1em',
            'y' : '1em',
            'inner sep': '0pt'}
    return tikz.get_tikz_response([], '', tikz_body, fs.tikzformat, options)

