"""Look at the eigendecomposition of the Laplacian of a decorated tree.

To each tip of the original tree are added N new leaves.
These new leaves have nonzero distance to the original tips of the tree.
The original tips of the tree are no longer tips
of the decorated tree.
"""


from StringIO import StringIO
import math

import numpy as np

from SnippetUtil import HandlingError
import SnippetUtil
import Form
import FormOut
import NewickIO
import Euclid
import FelTree
import const
import MatrixUtil

g_tree_string = const.read('20100730g').rstrip()


def get_form():
    """
    @return: a list of form objects
    """
    # define the form objects
    form_objects = [
            Form.MultiLine('tree_string', 'newick tree',
                g_tree_string),
            Form.Integer('N', 'tip duplication factor (N)', 3),
            Form.RadioGroup('weight_options', 'weight options', [
                Form.RadioItem('weight_n', 'N', True),
                Form.RadioItem('weight_sqrt_n', 'sqrt(N)')]),
            Form.CheckGroup('options', 'output options', [
                Form.CheckItem('lap', 'show laplacian', True),
                Form.CheckItem('eigvals', 'show eigenvalues', True),
                Form.CheckItem('eigvecs', 'show eigenvectors', True),
                Form.CheckItem('compare', 'compare leaf valuations', True)])]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):
    # get the tree
    tree = NewickIO.parse(fs.tree_string, FelTree.NewickTree)
    # get information about the tree topology
    internal = [id(node) for node in tree.gen_internal_nodes()]
    tips = [id(node) for node in tree.gen_tips()]
    vertices = internal + tips
    ntips = len(tips)
    ninternal = len(internal)
    nvertices = len(vertices)
    # get the ordered ids with the leaves first
    ordered_ids = vertices
    # get the full weighted adjacency matrix
    A = np.array(tree.get_affinity_matrix(ordered_ids))
    # compute the weighted adjacency matrix of the decorated tree
    p = ninternal
    q = ntips
    N = fs.N
    if fs.weight_n:
        weight = float(N)
    elif fs.weight_sqrt_n:
        weight = math.sqrt(N)
    A_aug = get_A_aug(A, weight, p, q, N)
    # compute the weighted Laplacian matrix of the decorated tree
    L_aug = Euclid.adjacency_to_laplacian(A_aug)
    # compute the eigendecomposition
    w, vt = np.linalg.eigh(L_aug)
    # show the output
    np.set_printoptions(linewidth=1000, threshold=10000)
    out = StringIO()
    if fs.lap:
        print >> out, 'Laplacian of the decorated tree:'
        print >> out, L_aug
        print >> out
    if fs.eigvals:
        print >> out, 'eigenvalues:'
        for x in w:
            print >> out, x
        print >> out
    if fs.eigvecs:
        print >> out, 'eigenvector matrix:'
        print >> out, vt
        print >> out
    if fs.compare:
        # get the distance matrix for only the original tips
        D_tips = np.array(tree.get_partial_distance_matrix(tips))
        X_tips = Euclid.edm_to_points(D_tips)
        # wring the approximate points out of the augmented tree
        X_approx = vt[p:p+q].T[1:1+q-1].T / np.sqrt(w[1:1+q-1])
        # do the comparison
        print >> out, 'points from tip-only MDS:'
        print >> out, X_tips
        print >> out
        print >> out, 'approximate points from decorated tree:'
        print >> out, X_approx
        print >> out
    return out.getvalue()

def get_A_aug(A, weight, p, q, N):
    """
    @param A: the weighted adjacency matrix for the full tree
    @param weight: reciprocal of lengths of decorating branches
    @param p: the number of internal nodes
    @param q: the number of tips
    @param N: the tip duplication factor
    @return: the weighted adjacency matrix for the decorated tree
    """
    k = p + q + N*q
    A_aug = np.zeros((k, k), dtype=float)
    # copy the original adjacency matrix
    for i in range(p+q):
        for j in range(p+q):
            A_aug[i,j] = A[i,j]
    # add the weighted connections between original tips and new leaves
    for a in range(N):
        for i in range(q):
            A_aug[p+q+q*a+i, p+i] += weight
            A_aug[p+i, p+q+q*a+i] += weight
    return A_aug
