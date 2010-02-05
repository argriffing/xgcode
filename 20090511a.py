"""Numerically verify that a certain function of a tree gives its full distance matrix.

The function of the tree is 2*(A + y*d*d' - deg)^-1.
Here A is the adjacency matrix of a full tree,
where edge weights are reciprocals of branch lengths.
The diagonal matrix deg consists of vertex degrees implied by A.
The constant y is the sum of branch lengths on the tree.
The column vector d and its transpose d' have elements -1 and 1 depending
on whether the corresponding vertex is an internal or a pendent vertex respectively.
If distance matrices are shown,
their rows and columns will be ordered consistently but arbitrarily.
"""

from StringIO import StringIO

import numpy

from SnippetUtil import HandlingError
import MatrixUtil
import Form
import NewickIO
import FelTree

def get_form():
    """
    @return: the body of a form
    """
    # define the default tree string
    tree_string = '(((a:0.05, b:0.05)ab:0.15, c:0.2)abc:0.8, x:1.0, (((m:0.05, n:0.05)mn:0.15, p:0.2)mnp:0.8, y:1.0)mnpy:1.0)abcxmnpy;'
    tree = NewickIO.parse(tree_string, FelTree.NewickTree)
    formatted_tree_string = NewickIO.get_narrow_newick_string(tree, 60)
    # define the form objects
    form_objects = [
            Form.MultiLine('tree', 'newick tree with branch lengths', formatted_tree_string),
            Form.CheckGroup('options', 'output options', [
                Form.CheckItem('show_direct_d', 'show the directly calculated full distance matrix', True),
                Form.CheckItem('show_clever_d', 'show the cleverly calculated full distance matrix', True),
                Form.CheckItem('show_closeness', 'show whether or not distance matrices are close', True)])]
    return form_objects


def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # get the tree
    tree = NewickIO.parse(fs.tree, FelTree.NewickTree)
    # get the arbitrarily ordered names
    ordered_names = set(node.get_name() for node in tree.preorder())
    # get the corresponding ordered ids
    name_to_id = dict((node.get_name(), id(node)) for node in tree.preorder())
    ordered_ids = [name_to_id[name] for name in ordered_names]
    # get the full distance matrix
    D_direct = numpy.array(tree.get_full_distance_matrix(ordered_ids))
    # get the full weighted adjacency matrix
    A = numpy.array(tree.get_affinity_matrix(ordered_ids))
    # get the full degree matrix
    degree_matrix = numpy.diag(numpy.sum(A, 0))
    # get the sum of the branch lengths
    n = len(ordered_names)
    gamma_inv = 0
    for i in range(n):
        for j in range(n):
            if i < j:
                if A[i][j]:
                    gamma_inv += 1.0 / A[i][j]
    gamma = 1.0 / gamma_inv
    # get the delta vector
    delta_list = []
    for row in A:
        nonzero_edge_count = sum(1 for x in row if x)
        delta_list.append(2 - nonzero_edge_count)
    d = numpy.array(delta_list)
    # get the full distance matrix using the clever formula
    J = numpy.ones((n, n))
    D_clever = 2*numpy.linalg.inv(A + gamma * numpy.outer(d, d) - degree_matrix)
    # check whether the distance matrices are close
    closeness_string = 'the distance matrices are close'
    if not numpy.allclose(D_direct, D_clever):
        closeness_string = 'the distance matrices are not close'
    # define the response
    out = StringIO()
    paragraphs = []
    if fs.show_direct_d:
        paragraph = [
                'directly calculated distance matrix:',
                MatrixUtil.m_to_string(D_direct)]
        paragraphs.append(paragraph)
    if fs.show_clever_d:
        paragraph = [
                'cleverly calculated distance matrix:',
                MatrixUtil.m_to_string(D_clever)]
        paragraphs.append(paragraph)
    if fs.show_closeness:
        paragraph = [
                'closeness:',
                closeness_string]
        paragraphs.append(paragraph)
    print >> out, '\n\n'.join('\n'.join(paragraph) for paragraph in paragraphs)
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
