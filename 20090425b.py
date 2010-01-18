"""Get a contrast matrix from a rooted tree.

Each column of the output matrix should be a contrast.
Elements in the output matrix with absolute values smaller than epsilon will be zeroed.
"""

import StringIO

from SnippetUtil import HandlingError
import Util
import NewickIO
import FelTree
import MatrixUtil
import Contrasts
import Form

def get_form():
    """
    @return: the body of a form
    """
    # define the default tree string
    tree_string = '((((a:0.05, b:0.05):0.15, c:0.2):0.8, x:1.0):0.5, (((m:0.05, n:0.05):0.15, p:0.2):0.8, y:1.0):0.5);'
    tree = NewickIO.parse(tree_string, FelTree.NewickTree)
    formatted_tree_string = NewickIO.get_narrow_newick_string(tree, 60)
    # define the ordered labels
    ordered_labels = list('abcxmnpy')
    # define the form objects
    form_objects = [
            Form.MultiLine('tree', 'rooted newick tree with branch lengths', formatted_tree_string),
            Form.MultiLine('labels', 'ordered labels', '\n'.join(ordered_labels)),
            Form.Float('epsilon', 'epsilon', '1e-10'),
            Form.RadioGroup('matrix_format', 'output matrix format', [
                Form.RadioItem('plain_format', 'plain', True),
                Form.RadioItem('r_format', 'R'),
                Form.RadioItem('matlab_format', 'MATLAB')])]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # get the tree
    tree = NewickIO.parse(fs.tree, FelTree.NewickTree)
    # read the ordered labels
    ordered_labels = list(Util.stripped_lines(StringIO.StringIO(fs.labels)))
    # validate the input
    observed_label_set = set(node.get_name() for node in tree.gen_tips())
    if set(ordered_labels) != observed_label_set:
        raise HandlingError('the ordered labels should match the labels of the leaves of the tree')
    # get the matrix of pairwise distances among the tips
    C = Contrasts.get_contrast_matrix(tree, ordered_labels)
    # set elements with small absolute value to zero
    C[abs(C) < fs.epsilon] = 0
    # start to prepare the reponse
    out = StringIO.StringIO()
    if fs.plain_format:
        print >> out, MatrixUtil.m_to_string(C)
    elif fs.matlab_format:
        print >> out, MatrixUtil.m_to_matlab_string(C)
    elif fs.r_format:
        print >> out, MatrixUtil.m_to_R_string(C)
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
