"""Follow the steps of a recursive tree reconstruction algorithm.

The exact bipartition criterion is a matrix function by Eric Stone.
"""

from StringIO import StringIO

import numpy as np

from SnippetUtil import HandlingError
import Util
import MatrixUtil
import NewickIO
import FelTree
import Clustering
import NeighborhoodJoining
from Form import RadioItem
import Form
import FormOut

def get_form():
    """
    @return: the body of a form
    """
    # define the default distance matrix and the ordered labels
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
    # define the default newick tree string
    tree_string = '(((a, b), c), x, (((m, n), p), y));'
    tree = NewickIO.parse(tree_string, FelTree.NewickTree)
    formatted_tree_string = NewickIO.get_narrow_newick_string(tree, 60)
    # define the form objects
    form_objects = [
            Form.Matrix('matrix', 'perturbed distance matrix',
                D, MatrixUtil.assert_predistance),
            Form.MultiLine('labels', 'ordered labels',
                '\n'.join(labels)),
            Form.MultiLine('tree', 'original tree (branch lengths optional)',
                formatted_tree_string),
            Form.RadioGroup('criterion', 'bipartition function', [
                RadioItem('exact', 'exact criterion', True),
                RadioItem('sign', 'spectral sign approximation'),
                RadioItem('threshold', 'spectral threshold approximation'),
                RadioItem('nj', 'neighbor joining criterion')]),
            Form.RadioGroup('recourse', 'recourse for degenerate partitions', [
                RadioItem('njrecourse', 'neighbor joining', True),
                RadioItem('halvingrecourse', 'leaf stem length halving')])]
    return form_objects

def get_form_out():
    return FormOut.Report()

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
    if len(set(ordered_labels)) != len(ordered_labels):
        raise HandlingError('the ordered labels must be unique')
    # read the criterion string, creating the splitter object
    if fs.exact:
        splitter = Clustering.StoneExactDMS()
    elif fs.sign:
        splitter = Clustering.StoneSpectralSignDMS()
    elif fs.threshold:
        splitter = Clustering.StoneSpectralThresholdDMS()
    elif fs.nj:
        splitter = Clustering.NeighborJoiningDMS()
    # Make sure that the splitter object
    # is appropriate for the size of the distance matrix.
    if splitter.get_complexity(len(D)) > 1000000:
        msg = 'use a smaller distance matrix or a faster bipartition function'
        raise HandlingError(msg)
    # read the original tree
    tree = NewickIO.parse(fs.tree, FelTree.NewickTree)
    if len(ordered_labels) != len(list(tree.gen_tips())):
        msg_a = 'the number of ordered labels should be the same '
        msg_b = 'as the number of tips in the tree'
        raise HandlingError(msg_a + msg_b)
    tree_tip_names = set(tip.name for tip in tree.gen_tips())
    if tree_tip_names != set(ordered_labels):
        msg_a = 'the leaf labels of the tree do not match '
        msg_b = 'the ordered labels of the distance matrix rows'
        raise HandlingError(msg_a + msg_b)
    # create the tree builder
    tree_builder = NeighborhoodJoining.ValidatingTreeBuilder(
            D.tolist(), ordered_labels, splitter)
    # Read the recourse string and set the corresponding method
    # in the tree builder.
    if fs.njrecourse:
        tree_builder.set_fallback_name('nj')
    elif fs.halvingrecourse:
        tree_builder.set_fallback_name('halving')
    # define the response
    out = StringIO()
    # set parameters of the tree validating tree builder
    tree_builder.set_original_tree(tree)
    tree_builder.set_output_stream(out)
    tree = tree_builder.build()
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
