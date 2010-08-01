"""Given a newick tree, calculate various matrices.
"""

from StringIO import StringIO

from scipy import linalg
import numpy as np

from SnippetUtil import HandlingError
import MatrixUtil
import Util
import NewickIO
import FelTree
import NeighborJoining
import Form
import FormOut
import iterutils

def get_form():
    """
    @return: the body of a form
    """
    # define the default tree string and ordered labels
    default_tree_string = '(a:1, (b:2, d:5):1, c:4);'
    tree = NewickIO.parse(default_tree_string, FelTree.NewickTree)
    ordered_labels = list(sorted(tip.name for tip in tree.gen_tips()))
    # define the form objects
    form_objects = [
            Form.MultiLine('tree', 'newick tree', default_tree_string),
            Form.MultiLine('inlabels', 'ordered labels',
                '\n'.join(ordered_labels)),
            Form.CheckGroup('options', 'output options', [
                Form.CheckItem('distance', 'path resistance matrix', True),
                Form.CheckItem('edge', 'edge resistance matrix'),
                Form.CheckItem('affinity', 'affinity matrix'),
                Form.CheckItem('laplacian', 'laplacian matrix'),
                Form.CheckItem('neglaplacian', 'negative laplacian matrix'),
                Form.CheckItem('Q', 'neighbor-joining Q matrix'),
                Form.CheckItem('labels', 'ordered labels', True)])]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # get the tree
    tree = NewickIO.parse(fs.tree, FelTree.NewickTree)
    # get the names of the tips of the tree
    alphabetically_ordered_states = list(sorted(node.name
        for node in tree.gen_tips()))
    n = len(alphabetically_ordered_states)
    if n < 2:
        raise HandlingError('the newick tree should have at least two leaves')
    # read the ordered labels
    states = []
    if fs.inlabels:
        states = Util.get_stripped_lines(fs.inlabels.splitlines())
    if len(states) > 1:
        if set(states) != set(alphabetically_ordered_states):
            msg_a = 'if ordered labels are provided, '
            msg_b = 'each should correspond to a leaf of the newick tree'
            raise HandlingError(msg_a + msg_b)
    else:
        states = alphabetically_ordered_states
    # create the distance matrix
    D = tree.get_distance_matrix(states)
    # create the laplacian matrix
    M = np.array(D)
    P = np.eye(n) - np.ones((n,n))/n
    L_pinv = - 0.5 * np.dot(P, np.dot(M, P))
    L = linalg.pinv(L_pinv)
    # start collecting the paragraphs
    paragraphs = []
    # show the distance matrix if requested
    if fs.distance:
        paragraph = StringIO()
        print >> paragraph, 'path resistance (distance) matrix:'
        print >> paragraph, MatrixUtil.m_to_string(D)
        paragraphs.append(paragraph.getvalue().strip())
    # show the edge matrix if requested
    if fs.edge:
        paragraph = StringIO()
        print >> paragraph, 'edge resistance matrix:'
        edge_matrix = L.copy()
        for i in range(n):
            for j in range(n):
                if i == j:
                    edge_matrix[i][j] = 0
                else:
                    edge_matrix[i][j] = -1.0 / edge_matrix[i][j]
        print >> paragraph, MatrixUtil.m_to_string(edge_matrix)
        paragraphs.append(paragraph.getvalue().strip())
    # show the affinity matrix if requested
    if fs.affinity:
        paragraph = StringIO()
        print >> paragraph, 'affinity matrix:'
        affinity_matrix = L.copy()
        for i in range(n):
            for j in range(n):
                if i == j:
                    affinity_matrix[i][j] = 0
                else:
                    affinity_matrix[i][j] *= -1
        print >> paragraph, MatrixUtil.m_to_string(affinity_matrix)
        paragraphs.append(paragraph.getvalue().strip())
    # show the laplacian matrix if requested
    if fs.laplacian:
        paragraph = StringIO()
        print >> paragraph, 'laplacian matrix:'
        print >> paragraph, MatrixUtil.m_to_string(L)
        paragraphs.append(paragraph.getvalue().strip())
    # show the negative laplacian matrix if requested
    if fs.neglaplacian:
        paragraph = StringIO()
        print >> paragraph, 'negative laplacian matrix:'
        print >> paragraph, MatrixUtil.m_to_string(-L)
        paragraphs.append(paragraph.getvalue().strip())
    # show the neighbor joining Q matrix
    if fs.Q:
        Q = NeighborJoining.get_Q_matrix(D)
        paragraph = StringIO()
        print >> paragraph, 'neighbor-joining Q matrix:'
        print >> paragraph, MatrixUtil.m_to_string(Q)
        paragraphs.append(paragraph.getvalue().strip())
    # show the ordered labels if requested
    if fs.labels:
        paragraph = StringIO()
        print >> paragraph, 'ordered labels:'
        print >> paragraph, '\n'.join(states)
        paragraphs.append(paragraph.getvalue().strip())
    # prepare the reponse
    text = '\n\n'.join(paragraphs) + '\n'
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, text
