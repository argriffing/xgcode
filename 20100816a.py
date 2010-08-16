"""Relate the Schur complement, tree tip duplication, and the Fiedler cut.

The idea is that duplicating the tips of a tree does not
really affect the scaled Fiedler vector of the Schur complement of
the internal nodes in the Laplacian.
"""

from StringIO import StringIO
import math

import numpy as np

from SnippetUtil import HandlingError
import Form
import FormOut
import NewickIO
import Euclid
import FelTree
import SchurAlgebra
import EigUtil
import const

g_tree_string = const.read('20100730g').rstrip()

def get_form():
    """
    @return: a list of form objects
    """
    # define the form objects
    form_objects = [
            Form.MultiLine('tree_string', 'newick tree',
                g_tree_string),
            Form.Float('strength', 'duplicate node pair connection strength',
                100)]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_augmented_adjacency(A, ntips, strength):
    """
    @param A: adjacency matrix
    @param ntips: the number of tips of the tree defining the adjacency matrix
    @param strength: the strength of the duplicate pair connection
    @return: a new adjacency matrix with tip extensions at the end
    """
    n = len(A)
    nnew = n + ntips
    ninternal = n - ntips
    P = np.zeros((nnew, nnew))
    for i in range(nnew):
        for j in range(nnew):
            if i < n and j < n:
                P[i,j] = A[i,j]
            elif i < ntips and j == n + i:
                P[i,j] = strength
            elif i == n + j and j < ntips:
                P[i,j] = strength
    return P

def get_response_content(fs):
    # build the newick tree from the string
    tree = NewickIO.parse(fs.tree_string, FelTree.NewickTree)
    nvertices = len(list(tree.preorder()))
    nleaves = len(list(tree.gen_tips()))
    # get ordered ids with the leaves first
    ordered_ids = get_ordered_ids(tree)
    # get the adjacency matrix and the augmented adjacency matrix
    A = np.array(tree.get_affinity_matrix(ordered_ids))
    A_aug = get_augmented_adjacency(A, nleaves, fs.strength)
    # get the laplacian matrices
    L = Euclid.adjacency_to_laplacian(A)
    L_aug = Euclid.adjacency_to_laplacian(A_aug)
    # get the schur complements
    R = SchurAlgebra.mschur(L, set(range(nleaves, nvertices)))
    R_aug = SchurAlgebra.mschur(L_aug, set(range(nleaves, nvertices)))
    # get the scaled Fiedler vectors
    w, v = EigUtil.principal_eigh(np.linalg.pinv(R))
    fiedler = v * math.sqrt(w)
    w_aug, v_aug = EigUtil.principal_eigh(np.linalg.pinv(R_aug))
    fiedler_aug = v_aug * math.sqrt(w_aug)
    # report the results
    np.set_printoptions(linewidth=200)
    out = StringIO()
    print >> out, 'Laplacian matrix:'
    print >> out, L
    print >> out
    print >> out, 'Schur complement of Laplacian matrix:'
    print >> out, R
    print >> out
    print >> out, 'scaled Fiedler vector:'
    print >> out, fiedler
    print >> out
    print >> out, 'augmented Laplacian matrix:'
    print >> out, L_aug
    print >> out
    print >> out, 'Schur complement of augmented Laplacian matrix:'
    print >> out, R_aug
    print >> out
    print >> out, 'scaled Fiedler vector of augmented matrix:'
    print >> out, fiedler_aug
    print >> out
    return out.getvalue()

def get_ordered_ids(tree):
    """
    Maybe I could use postorder here instead.
    @param tree: a tree
    @return: a list of ids beginning with the leaves
    """
    ordered_ids = []
    ordered_ids.extend(id(node) for node in tree.gen_tips())
    ordered_ids.extend(id(node) for node in tree.gen_internal_nodes())
    return ordered_ids
