"""Given a tree, use the normalized Laplacian to define a node distribution.

Given a newick tree, use the normalized Laplacian
to define a distribution over the nodes.
"""

from StringIO import StringIO
import math

import numpy as np

from SnippetUtil import HandlingError
import Form
import FormOut
import MatrixUtil
import NewickIO
import FelTree
import HtmlTable
import const

g_default_string = const.read('20100730m')

def get_form():
    """
    @return: a list of form objects
    """
    # define the default tree string
    tree = NewickIO.parse(g_default_string, FelTree.NewickTree)
    formatted_tree_string = NewickIO.get_narrow_newick_string(tree, 60)
    # return the form objects
    return [
            Form.MultiLine('tree', 'newick tree with branch lengths',
                formatted_tree_string),
            Form.Integer('precision', 'precision', 4, low=2, high=17)]

def get_form_out():
    return FormOut.Html()

def get_response(fs):
    """
    @param fs: a FieldStorage object decorated with field values
    @return: a (response_headers, response_text) pair
    """
    # get the tree
    tree = NewickIO.parse(fs.tree, FelTree.NewickTree)
    # define the ordering of the nodes
    ordered_names = list(sorted(node.get_name() for node in tree.preorder()))
    name_to_id = dict((node.get_name(), id(node)) for node in tree.preorder())
    ordered_ids = [name_to_id[name] for name in ordered_names]
    # get the affinity matrix
    A = np.array(tree.get_affinity_matrix(ordered_ids))
    # get the normalized laplacian
    for row, name in zip(A, ordered_names):
        assert sum(row), name
    row_sums = [sum(row) for row in A]
    n = len(A)
    L_script = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            if i == j:
                L_script[i][j] = 1
            L_script[i][j] -= A[i][j] / math.sqrt(row_sums[i]*row_sums[j])
    # get the eigendecomposition
    eigenvalues, eigenvector_transposes = np.linalg.eigh(L_script)
    eigenvectors = eigenvector_transposes.T
    eigensystem = [(abs(w), w, v.tolist()) for w, v in zip(eigenvalues, eigenvectors)]
    sorted_eigensystem = list(reversed(sorted(eigensystem)))
    sorted_abs_eigenvalues, sorted_eigenvalues, sorted_eigenvectors = zip(*sorted_eigensystem)
    M = zip(*sorted_eigenvectors)
    # turn the numbers into strings
    format_string = '%%.%df' % fs.precision
    M_strings = [[format_string % value for value in row] for row in M]
    sorted_eigenvalue_strings = [format_string % w for w in sorted_eigenvalues]
    # write the html
    out = StringIO()
    print >> out, '<html>'
    print >> out, '<body>'
    print >> out, HtmlTable.get_labeled_table_string(sorted_eigenvalue_strings, ordered_names, M_strings)
    print >> out, '</body>'
    print >> out, '</html>'
    # write the response
    response_headers = [('Content-Type', 'text/html')]
    return response_headers, out.getvalue().strip()
