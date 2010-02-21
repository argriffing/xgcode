"""Make a newick tree from a distance matrix by splitting instead of joining.
"""

from StringIO import StringIO

import numpy as np

from SnippetUtil import HandlingError
import Util
import MatrixUtil
import Clustering
import NeighborhoodJoining
import NewickIO
import Form

def get_form():
    """
    @return: the body of a form
    """
    # define the default distance matrix
    # this is from figure two of a paper called why neighbor joining works
    D = np.array([
        [  0, 2.7, 2.6, 2.6, 2.6, 4.4, 4.4, 4.4],
        [2.7,   0, 4.4, 4.4, 4.4, 2.6, 2.6, 2.6],
        [2.6, 4.4,   0, 0.1, 0.4, 2.7, 2.7, 2.7],
        [2.6, 4.4, 0.1,   0, 0.4, 2.7, 2.7, 2.7],
        [2.6, 4.4, 0.4, 0.4,   0, 2.7, 2.7, 2.7],
        [4.4, 2.6, 2.7, 2.7, 2.7,   0, 0.1, 0.4],
        [4.4, 2.6, 2.7, 2.7, 2.7, 0.1,   0, 0.4],
        [4.4, 2.6, 2.7, 2.7, 2.7, 0.4, 0.4,   0]])
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
    ordered_labels = Util.get_stripped_lines(StringIO(fs.labels))
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
    out = StringIO()
    print >> out, NewickIO.get_newick_string(tree)
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
