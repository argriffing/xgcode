"""Use neighbor joining to make a newick tree from a distance matrix.
"""

from StringIO import StringIO

import numpy as np

from SnippetUtil import HandlingError
import Util
import MatrixUtil
import NewickIO
import NeighborJoining
import Form
import FormOut

def get_form():
    """
    @return: the body of a form
    """
    # define the default distance matrix
    # this is from figure 2 of a paper called why neighbor joining works
    D = np.array([
        [0.0, 3.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0],
        [3.0, 0.0, 3.0, 3.0, 3.0, 2.0, 2.0, 2.0],
        [2.0, 3.0, 0.0, 0.1, 0.4, 3.0, 3.0, 3.0],
        [2.0, 3.0, 0.1, 0.0, 0.4, 3.0, 3.0, 3.0],
        [2.0, 3.0, 0.4, 0.4, 0.0, 3.0, 3.0, 3.0],
        [3.0, 2.0, 3.0, 3.0, 3.0, 0.0, 0.1, 0.4],
        [3.0, 2.0, 3.0, 3.0, 3.0, 0.1, 0.0, 0.4],
        [3.0, 2.0, 3.0, 3.0, 3.0, 0.4, 0.4, 0.0]])
    labels = list('xyabcmnp')
    # define the form objects
    form_objects = [
            Form.Matrix('matrix', 'distance matrix',
                D, MatrixUtil.assert_predistance),
            Form.MultiLine('labels', 'ordered labels',
                '\n'.join(labels))]
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
    ordered_labels = Util.get_stripped_lines(StringIO(fs.labels))
    if len(ordered_labels) != len(D):
        msg_a = 'the number of ordered labels should be the same '
        msg_b = 'as the number of rows in the matrix'
        raise HandlingError(msg_a + msg_b)
    # get the newick tree
    tree = NeighborJoining.make_tree(D.tolist(), ordered_labels)
    # define the response
    out = StringIO()
    print >> out, NewickIO.get_newick_string(tree)
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
