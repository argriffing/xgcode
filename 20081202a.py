"""Given a newick tree, show the eigendecomposition of the Gower matrix.

Here the Gower matrix is defined as -(1/2)HDH where D is the
distance matrix implied by the tree and H is the conformant centering matrix.
The column labels are the eigenvalues,
and the other numbers are the eigenvectors.
"""

from StringIO import StringIO

import numpy as np

from SnippetUtil import HandlingError
import Form
import FormOut
import MatrixUtil
import NewickIO
import FelTree
import HtmlTable
import const

g_default_string = const.read('20100730q')

def get_form():
    """
    @return: a list of form objects
    """
    # define the default tree string
    tree = NewickIO.parse(g_default_string, FelTree.NewickTree)
    formatted_tree_string = NewickIO.get_narrow_newick_string(tree, 60)
    # return the form objects
    form_objects = [
            Form.MultiLine('tree', 'newick tree with branch lengths',
                formatted_tree_string)]
    return form_objects

def get_form_out():
    return FormOut.Html()

def get_response(fs):
    """
    @param fs: a FieldStorage object decorated with field values
    @return: a (response_headers, response_text) pair
    """
    # get the tree
    tree = NewickIO.parse(fs.tree, FelTree.NewickTree)
    ordered_names = list(sorted(node.name for node in tree.gen_tips()))
    n = len(ordered_names)
    if n < 2:
        raise HandlingError('the newick tree should have at least two leaves')
    # get the eigendecomposition
    D = np.array(tree.get_distance_matrix(ordered_names))
    G = (-0.5) * MatrixUtil.double_centered(D)
    eigenvalues, eigenvector_transposes = np.linalg.eigh(G)
    eigenvectors = eigenvector_transposes.T
    sorted_eigensystem = list(reversed(list(sorted((w, v) for w, v in zip(eigenvalues, eigenvectors)))))
    sorted_eigenvalues, sorted_eigenvectors = zip(*sorted_eigensystem)
    M = zip(*sorted_eigenvectors)
    # write the html
    out = StringIO()
    print >> out, '<html>'
    print >> out, '<body>'
    print >> out, HtmlTable.get_labeled_table_string(sorted_eigenvalues, ordered_names, M)
    print >> out, '</body>'
    print >> out, '</html>'
    # write the response
    response_headers = [('Content-Type', 'text/html')]
    return response_headers, out.getvalue().strip()
