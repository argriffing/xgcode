"""Visualize an internal state of a tree reconstruction algorithm.

Visualize an internal state of a distance based tree reconstruction algorithm.
"""

from StringIO import StringIO
import tempfile

import numpy as np

from SnippetUtil import HandlingError
import SnippetUtil
import Form
import Util
import MatrixUtil
import GraphNHJ
import NeighborJoining

def get_form():
    """
    @return: the body of a form
    """
    # define the default distance matrix
    D = np.array(NeighborJoining.g_mito_matrix)
    ordered_labels = ['gorilla', 'orangutan', 'human', 'chimp', 'gibbon']
    # get some strings that are default values
    ordered_label_string = '\n'.join(ordered_labels)
    # create some form objects
    form_objects = [
            Form.Matrix('matrix', 'distance matrix',
                D, MatrixUtil.assert_predistance),
            Form.MultiLine('labels', 'ordered labels',
                ordered_label_string),
            Form.Integer('iteration', 'the iteration to visualize',
                1),
            Form.RadioGroup('contentdisposition', 'image display options', [
                Form.RadioItem('inline', 'view the image', True),
                Form.RadioItem('attachment', 'download the image')])]
    # return the objects
    return form_objects

def get_image_string(D, ordered_labels, iteration):
    """
    @param D: a numpy distance matrix
    @param ordered_labels: each label is the name of a distance matrix row
    @param iteration: the iteration whose graphical view is requested
    @return: the string representation of a png graphics file
    """
    # do some basic validation of the input
    min_iteration = 1
    max_iteration = len(D) - 2
    if not (min_iteration <= iteration <= max_iteration):
        raise ValueError('invalid iteration')
    # convert the iteration to the number of iterations to do
    niterations = iteration - 1
    # create the initial graph corresponding to the distance matrix
    treegraph = GraphNHJ.Tree(D.tolist())
    # do some iterations of graph transformation
    for i in range(niterations):
        node = treegraph.get_decomposible_node()
        if node:
            treegraph.decompose(node, GraphNHJ.get_augmented_gower_selection)
        else:
            ValueError('no decomposible node was found')
    # define the temporary filename to be used by pygraphviz
    filename = tempfile.mktemp() + '.png'
    # get the pygraphviz graph that we want to draw
    G = treegraph.get_pygraphviz_graph(ordered_labels)
    # try to prevent the nodes from overlapping when they are drawn
    G.graph_attr['overlap'] = 'false'
    # draw the figure
    G.layout()
    G.draw(filename)
    # read the binary image string from the image file
    fin = open(filename, 'rb')
    image_string = fin.read()
    fin.close()
    # return the image string
    return image_string

def get_response(fs):
    """
    @param fs: a FieldStorage object decorated with field values
    @return: a (response_headers, response_text) pair
    """
    # read the matrix
    D = fs.matrix
    if len(D) < 3:
        msg = 'the distance matrix should have at least three rows'
        raise HandlingError(msg)
    # read the ordered labels
    ordered_labels = Util.get_stripped_lines(StringIO(fs.labels))
    if not ordered_labels:
        raise HandlingError('no ordered labels were provided')
    if len(ordered_labels) != len(D):
        msg_a = 'the number of ordered labels '
        msg_b = 'should be the same as the number of rows in the matrix'
        raise HandlingError(msg_a + msg_b)
    if len(set(ordered_labels)) != len(ordered_labels):
        raise HandlingError('the ordered labels must be unique')
    # read the index of the iteration that will be visualized
    min_iteration = 1
    max_iteration = len(D) - 2
    iteration = fs.iteration
    if not (min_iteration <= iteration <= max_iteration):
        msg_a = 'the iteration index '
        msg_b = 'should be in [%d, %d]' % (min_iteration, max_iteration)
        raise HandlingError(msg_a + msg_b)
    # get the image string
    image_string = get_image_string(D, ordered_labels, iteration)
    # specify the content type (hard coded to png)
    response_headers = []
    response_headers.append(('Content-Type', 'image/png'))
    # specify the content disposition (specified by the user)
    response_headers.append(('Content-Disposition', "%s; filename=%s" % (fs.contentdisposition, 'graph.png')))
    # return the response
    return response_headers, image_string

