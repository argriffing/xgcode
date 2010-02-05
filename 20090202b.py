"""Given a newick tree, use the normalized Laplacian to define a distribution over the nodes.
"""

from StringIO import StringIO
import math

import numpy
from scipy import linalg

from SnippetUtil import HandlingError
import Form
import MatrixUtil
import NewickIO
import FelTree
import HtmlTable

def get_form():
    """
    @return: a list of form objects
    """
    # define the default tree string
    tree_string = '(((a:0.05, b:0.05)ab:0.15, c:0.2)abc:0.8, x:1.0, (((m:0.05, n:0.05)mn:0.15, p:0.2)mnp:0.8, y:1.0)mnpy:1.0)abcxmnpy;'
    tree = NewickIO.parse(tree_string, FelTree.NewickTree)
    formatted_tree_string = NewickIO.get_narrow_newick_string(tree, 60)
    # return the form objects
    return [
            Form.MultiLine('tree', 'newick tree with branch lengths', formatted_tree_string),
            Form.Integer('precision', 'precision', 4, low=2, high=17)]

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
    A = numpy.array(tree.get_affinity_matrix(ordered_ids))
    # get the normalized laplacian
    for row, name in zip(A, ordered_names):
        assert sum(row), name
    row_sums = [sum(row) for row in A]
    n = len(A)
    L_script = numpy.zeros((n,n))
    for i in range(n):
        for j in range(n):
            if i == j:
                L_script[i][j] = 1
            L_script[i][j] -= A[i][j] / math.sqrt(row_sums[i]*row_sums[j])
    # get the eigendecomposition
    eigenvalues, eigenvector_transposes = linalg.eigh(L_script)
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
