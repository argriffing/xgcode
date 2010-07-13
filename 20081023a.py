"""Visualize the steps of a tree inference algorithm. [UNFINISHED]

Visualize the steps of a distance based
tree reconstruction algorithm.
Output is a tarfile of a sequence of tree images.
"""

from StringIO import StringIO

import numpy as np

from SnippetUtil import HandlingError
import Util
import MatrixUtil
import NewickIO
import Clustering
import NeighborJoining
import NeighborhoodJoining
from Form import RadioItem
import Form
import FormOut

def get_form():
    """
    @return: the body of a form
    """
    # define the default distance matrix
    D = np.array(NeighborJoining.g_mito_matrix)
    # define the default label order
    labels = ['gorilla', 'orangutan', 'human', 'chimp', 'gibbon']
    # define the form objects
    form_objects = [
            Form.Matrix('matrix', 'distance matrix',
                D, MatrixUtil.assert_predistance),
            Form.MultiLine('labels', 'ordered labels',
                '\n'.join(labels)),
            Form.RadioGroup('criterion', 'tree reconstruction criterion', [
                RadioItem('nj_specific', 'neighbor joining (specific)'),
                RadioItem('nj_general', 'neighbor joining (general)'),
                RadioItem('sign', 'spectral sign approximation', True),
                RadioItem('random', 'random bipartition')])]
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
    if fs.sign:
        splitter = Clustering.StoneSpectralSignDMS()
    elif fs.threshold:
        splitter = Clustering.StoneSpectralThresholdDMS()
    elif fs.nj_general:
        splitter = Clustering.NeighborJoiningDMS()
    elif fs.nj_specific:
        splitter = None
    # Make sure that the splitter object is appropriate
    # for the size of the distance matrix.
    if splitter.get_complexity(len(D)) > 1000000:
        msg_a = 'use a smaller distance matrix '
        msg_b = 'or a faster bipartition function'
        raise HandlingError(msg_a + msg_b)
    # create the tree builder
    tree_builder = NeighborhoodJoining.TreeBuilder(
            D.tolist(), ordered_labels, splitter)
    tree_builder.set_fallback_name('nj')
    # define the response
    out = StringIO()
    # build the tree
    tree = tree_builder.build()
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
