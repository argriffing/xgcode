"""Get a contrast matrix from a rooted tree.

Each column of the output matrix should be a contrast.
Elements in the output matrix with absolute values smaller than epsilon
will be zeroed.
"""

from SnippetUtil import HandlingError
import Util
import NewickIO
import FelTree
import MatrixUtil
import Contrasts
import Form
import FormOut
import const

g_default_string = const.read('20100730n')

#FIXME matrix output formats

def get_form():
    """
    @return: the body of a form
    """
    # define the default tree string
    tree = NewickIO.parse(g_default_string, FelTree.NewickTree)
    formatted_tree_string = NewickIO.get_narrow_newick_string(tree, 60)
    # define the ordered labels
    ordered_labels = list('abcxmnpy')
    # define the form objects
    form_objects = [
            Form.MultiLine('tree', 'rooted newick tree with branch lengths',
                formatted_tree_string),
            Form.MultiLine('labels', 'ordered labels',
                '\n'.join(ordered_labels)),
            Form.Float('epsilon', 'epsilon', '1e-10'),
            Form.RadioGroup('matrix_format', 'output matrix format', [
                Form.RadioItem('plain_format', 'plain', True),
                Form.RadioItem('r_format', 'R'),
                Form.RadioItem('matlab_format', 'MATLAB')])]
    return form_objects

def get_form_out():
    return FormOut.ContextDependent()

def get_response_content(fs):
    # get the tree
    tree = NewickIO.parse(fs.tree, FelTree.NewickTree)
    # read the ordered labels
    ordered_labels = Util.get_stripped_lines(fs.labels.splitlines())
    # validate the input
    observed_label_set = set(node.get_name() for node in tree.gen_tips())
    if set(ordered_labels) != observed_label_set:
        msg = 'the labels should match the labels of the leaves of the tree'
        raise HandlingError(msg)
    # get the matrix of pairwise distances among the tips
    C = Contrasts.get_contrast_matrix(tree, ordered_labels)
    # set elements with small absolute value to zero
    C[abs(C) < fs.epsilon] = 0
    # return the reponse
    if fs.plain_format:
        return MatrixUtil.m_to_string(C) + '\n'
    elif fs.matlab_format:
        return MatrixUtil.m_to_matlab_string(C) + '\n'
    elif fs.r_format:
        return MatrixUtil.m_to_R_string(C) + '\n'
