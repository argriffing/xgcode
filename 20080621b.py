"""Make a newick tree from a distance matrix by splitting instead of joining.
"""

import numpy as np

from SnippetUtil import HandlingError
import Util
import MatrixUtil
import Clustering
import NeighborhoodJoining
import NewickIO
import Form
import FormOut
import const

# this is a perturbed distance matrix not an exact distance matrix
g_d = const.read('20100730p')

def get_form():
    """
    @return: the body of a form
    """
    lines = Util.get_stripped_lines(g_d.splitlines())
    D = np.array(MatrixUtil.read_matrix(lines))
    labels = list('xyabcmnp')
    # define the form objects
    form_objects = [
            Form.Matrix('matrix', 'distance matrix',
                D, MatrixUtil.assert_predistance),
            Form.MultiLine('labels', 'ordered labels',
                '\n'.join(labels)),
            Form.RadioGroup('recourse', 'recourse for degenerate partitions', [
                Form.RadioItem('njrecourse', 'neighbor joining', True),
                Form.RadioItem('halvingrecourse', 'stem length halving')])]
    return form_objects

def get_form_out():
    return FormOut.Newick()

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # read the matrix
    D = fs.matrix
    if len(D) < 3:
        raise HandlingError('the matrix should have at least three rows')
    # read the ordered labels
    ordered_labels = Util.get_stripped_lines(fs.labels.splitlines())
    if len(ordered_labels) != len(D):
        msg_a = 'the number of ordered labels should be the same '
        msg_b = 'as the number of rows in the matrix'
        raise HandlingError(msg_a + msg_b)
    # create the tree building object
    splitter = Clustering.StoneExactDMS()
    tree_builder = NeighborhoodJoining.TreeBuilder(
            D.tolist(), ordered_labels, splitter)
    # Read the recourse string and set the corresponding method
    # in the tree builder.
    recourse_string = fs.getfirst('recourse')
    if fs.njrecourse:
        tree_builder.set_fallback_name('nj')
    elif fs.halvingrecourse:
        tree_builder.set_fallback_name('halving')
    # assert that the computation will not take too long
    if tree_builder.get_complexity() > 1000000:
        raise HandlingError('this computation would take too long')
    # build the tree
    tree = tree_builder.build()
    # define the response
    text = NewickIO.get_newick_string(tree) + '\n'
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, text
