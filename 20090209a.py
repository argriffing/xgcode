"""Check a determinantal property of the full Laplacian matrix of a tree.

Note that this determinantal property fails for the full Laplacian matrix.
The determinant should be near zero
if the two sides define a valid split of the tree.
That is, the reported determinant should be near zero
when the reported branch length is positive,
and the reported determinant should be nonzero
when the reported branch length is negative.
"""

from StringIO import StringIO

import numpy as np

from SnippetUtil import HandlingError
import MatrixUtil
import Form
import FormOut
import NewickIO
import FelTree
import Clustering
import const

g_default_string = const.read('20100730m')

def get_form():
    """
    @return: the body of a form
    """
    # define the default tree string
    tree = NewickIO.parse(g_default_string, FelTree.NewickTree)
    formatted_tree_string = NewickIO.get_narrow_newick_string(tree, 60)
    # define the form objects
    form_objects = [
            Form.MultiLine('tree',
                'newick tree with branch lengths', formatted_tree_string),
            Form.SingleLine('lhs_a',
                'the first taxon on one side of the split', 'a'),
            Form.SingleLine('lhs_b',
                'the second taxon on one side of the split', 'b'),
            Form.SingleLine('rhs_a',
                'the first taxon on the other side of the split', 'x'),
            Form.SingleLine('rhs_b',
                'the second taxon on the other side of the split', 'y'),
            Form.CheckGroup('options', 'output options', [
                Form.CheckItem('show_response',
                    'show the full Laplacian matrix'),
                Form.CheckItem('show_reduced_response',
                    'show the 2x2 submatrix'),
                Form.CheckItem('show_blen',
                    'show the branch length implied by the split')])]
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
    # assert the the given labels are tips of the tree
    available_name_set = set(node.get_name() for node in tree.preorder())
    user_name_set = set([fs.lhs_a, fs.lhs_b, fs.rhs_a, fs.rhs_b])
    bad_names = user_name_set - available_name_set
    if bad_names:
        msg = 'these labels are not valid: %s' % ', '.join(bad_names)
        raise HandlingError(msg)
    # get the ordered names and ordered ids
    ordered_names = list(sorted(node.get_name() for node in tree.preorder()))
    name_to_id = dict((node.get_name(), id(node)) for node in tree.preorder())
    ordered_ids = [name_to_id[name] for name in ordered_names]
    # get the full affinity matrix
    A = tree.get_affinity_matrix(ordered_ids)
    # get the full Laplacian matrix
    row_sums = [sum(row) for row in A]
    n = len(ordered_names)
    L = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                L[i][i] = row_sums[i]
            else:
                L[i][j] = -A[i][j]
    # get the full distance matrix
    D = tree.get_full_distance_matrix(ordered_ids)
    # get the two by two matrix
    name_to_index = dict((name, i) for i, name in enumerate(ordered_names))
    L_reduced = np.zeros((2,2))
    la = name_to_index[fs.lhs_a]
    lb = name_to_index[fs.lhs_b]
    ra = name_to_index[fs.rhs_a]
    rb = name_to_index[fs.rhs_b]
    L_reduced[0][0] = L[la][ra]
    L_reduced[0][1] = L[la][rb]
    L_reduced[1][0] = L[lb][ra]
    L_reduced[1][1] = L[lb][rb]
    epsilon = 1e-13
    epsilon = 0.0000000000001
    criterion = np.linalg.det(L_reduced)
    if abs(criterion) < epsilon:
        criterion = 0
    # in analogy to the four point condition, use two different ways of calculating the distance
    blen_a = (D[la][rb] + D[lb][ra] - D[la][lb] - D[ra][rb]) / 2.0
    blen_b = (D[la][ra] + D[lb][rb] - D[la][lb] - D[ra][rb]) / 2.0
    blen = min(blen_a, blen_b)
    # define the response
    out = StringIO()
    paragraphs = []
    if fs.show_response:
        paragraph = [
                'full Laplacian matrix with rows ordered alphabetically by leaf label:',
                MatrixUtil.m_to_string(L)]
        paragraphs.append(paragraph)
    if fs.show_reduced_response:
        paragraph = [
                '2x2 submatrix of the Laplacian matrix:',
                MatrixUtil.m_to_string(L_reduced)]
        paragraphs.append(paragraph)
    if True:
        paragraph = [
                'determinant of the 2x2 submatrix of the Laplacian matrix:',
                str(criterion)]
        paragraphs.append(paragraph)
    if fs.show_blen:
        paragraph = [
                'branch length defined by the split:',
                str(blen)]
        paragraphs.append(paragraph)
    print >> out, '\n\n'.join('\n'.join(paragraph) for paragraph in paragraphs)
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
