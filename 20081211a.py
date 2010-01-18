"""Given a newick tree, compute a test matrix in three different ways.

Output matrix rows are ordered alphabetically by the leaf names.
This snippet is not so useful anymore,
but it helped find a mistake in some of my equations.
"""

import StringIO

from scipy import linalg
import numpy

from SnippetUtil import HandlingError
import Form
import MatrixUtil
import NewickIO
import FelTree
import Euclid
import Clustering

def get_form():
    """
    @return: a list of form objects
    """
    tree_string = '(a:2, b:2, (c:2, d:2):2);'
    form_item = Form.MultiLine('tree', 'newick tree with branch lengths', tree_string)
    return [form_item]

def transform_a(L, D_partial):
    """
    @param L: the full tree laplacian partitioned with the leaves before the internal nodes
    @param D_partial: distances between leaves
    @return: the S matrix
    """
    nleaves = len(D_partial)
    # get the Moore-Penrose inverse of the laplacian
    L_pinv = linalg.pinv(numpy.array(L))
    # extract the principal k by k matrix from L_pinv
    reduced_L_pinv = numpy.zeros((nleaves, nleaves))
    for i in range(nleaves):
        for j in range(nleaves):
            reduced_L_pinv[i][j] = L_pinv[i][j]
    # get the Moore-Penrose inverse of this reduced matrix
    S = linalg.pinv(reduced_L_pinv)
    return S

def transform_b(L, D_partial):
    """
    @param L: the full tree laplacian partitioned with the leaves before the internal nodes
    @param D_partial: distances between leaves
    @return: the S matrix
    """
    n = len(L)
    nleaves = len(D_partial)
    nancestors = len(L) - len(D_partial)
    # decompose the full laplacian
    L_out = numpy.zeros((nleaves, nleaves))
    for i in range(nleaves):
        for j in range(nleaves):
            L_out[i][j] = L[i][j]
    L_in = numpy.zeros((nancestors, nancestors))
    for i in range(nancestors):
        for j in range(nancestors):
            L_in[i][j] = L[nleaves+i][nleaves+j]
    A = numpy.zeros((nleaves, nancestors))
    for i in range(nleaves):
        for j in range(nancestors):
            A[i][j] = L[i][nleaves+j]
    # get the transformed matrix
    L_in_inv = linalg.inv(L_in)
    S = L_out - numpy.dot(A, numpy.dot(L_in_inv, A.T))
    return S

def transform_c(L, D_partial):
    """
    @param L: the full tree laplacian partitioned with the leaves before the internal nodes
    @param D_partial: distances between leaves
    @return: the S matrix
    """
    M = Clustering.get_R_stone(D_partial)
    S = -2*M
    return S

def get_response(fs):
    """
    @param fs: a decorated FieldStorage object
    @return: a (response_headers, response_text) pair
    """
    out = StringIO.StringIO()
    # get the tree
    tree = NewickIO.parse(fs.tree, FelTree.NewickTree)
    n = len(tree.preorder())
    # get the ids ordered so that the leaf ids are first
    leaf_name_id_pairs = [(node.get_name(), id(node)) for node in tree.gen_tips()]
    ordered_leaf_ids = [node_id for name, node_id in sorted(leaf_name_id_pairs)]
    all_ids = [id(node) for node in tree.preorder()]
    ordered_non_leaf_ids = list(set(all_ids) - set(ordered_leaf_ids))
    ordered_ids = ordered_leaf_ids + ordered_non_leaf_ids
    # get the full laplacian matrix with small values set to zero
    D_full = tree.get_full_distance_matrix(ordered_ids)
    L_full = Euclid.edm_to_laplacian(numpy.array(D_full))
    epsilon = 1e-10
    L_full[numpy.abs(L_full) < epsilon] = 0
    # get the partial distance matrix
    D_partial = tree.get_partial_distance_matrix(ordered_leaf_ids)
    # show the partial distance matrix
    print >> out, 'partial distance matrix (leaves only):'
    print >> out, MatrixUtil.m_to_string(D_partial)
    print >> out
    # show the full laplacian matrix
    print >> out, 'full laplacian matrix (leaves first):'
    print >> out, MatrixUtil.m_to_string(L_full)
    print >> out
    # show the output matrices
    functions = (transform_a, transform_b, transform_c)
    names = ('first transformation', 'second transformation', 'third transformation')
    for name, fn in zip(names, functions):
        S = fn(L_full, D_partial)
        print >> out, name + ':'
        print >> out, MatrixUtil.m_to_string(S)
        print >> out
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
