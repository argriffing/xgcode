"""Look at the eigendecomposition of a centered augmented distance matrix.

Tips of the tree are duplicated.
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
            Form.Integer('ndups', 'add this many extra duplicates per tip', 3)]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_augmented_distance(D, ntips, ndups):
    """
    @param D: distance matrix
    @param ntips: the number of tips of the original tree
    @param ndups: add this many extra duplicates per tip
    @return: a new distance matrix with tip extensions at the end
    """
    n = len(D)
    nnew = n + ntips*ndups
    ninternal = n - ntips
    P = np.zeros((nnew, nnew))
    for i in range(nnew):
        for j in range(nnew):
            if i < n and j < n:
                P[i,j] = D[i,j]
            elif i < n and j >= n:
                P[i,j] = D[i,(j-n)%ntips]
            elif j < n and i >= n:
                P[i,j] = D[(i-n)%ntips,j]
            else:
                P[i,j] = D[(i-n)%ntips,(j-n)%ntips]
    return P

def get_response_content(fs):
    # build the newick tree from the string
    tree = NewickIO.parse(fs.tree_string, FelTree.NewickTree)
    nvertices = len(list(tree.preorder()))
    nleaves = len(list(tree.gen_tips()))
    # get ordered ids with the leaves first
    ordered_ids = get_ordered_ids(tree)
    # get the distance matrix and the augmented distance matrix
    D = np.array(tree.get_partial_distance_matrix(ordered_ids))
    D_aug = get_augmented_distance(D, nleaves, fs.ndups)
    # get the laplacian matrix
    L = Euclid.edm_to_laplacian(D)
    # get the schur complement
    R = SchurAlgebra.mschur(L, set(range(nleaves, nvertices)))
    R_pinv = np.linalg.pinv(R)
    vals, vecs = EigUtil.eigh(R_pinv)
    # get the scaled Fiedler vector for the Schur complement
    w, v = EigUtil.principal_eigh(R_pinv)
    fiedler = v * math.sqrt(w)
    # get the eigendecomposition of the centered augmented distance matrix
    L_aug_pinv = Euclid.edm_to_dccov(D_aug)
    vals_aug, vecs_aug = EigUtil.eigh(L_aug_pinv)
    # get the scaled Fiedler vector for the augmented Laplacian
    w_aug, v_aug = EigUtil.principal_eigh(L_aug_pinv)
    fiedler_aug = v_aug * math.sqrt(w_aug)
    # report the results
    np.set_printoptions(linewidth=300)
    out = StringIO()
    print >> out, 'Laplacian matrix:'
    print >> out, L
    print >> out
    print >> out, 'Schur complement of Laplacian matrix:'
    print >> out, R
    print >> out
    print >> out, 'scaled Fiedler vector of Schur complement:'
    print >> out, fiedler
    print >> out
    print >> out, 'eigenvalues of pinv of Schur complement:'
    print >> out, vals
    print >> out
    print >> out, 'corresponding eigenvectors of pinv of Schur complement:'
    print >> out, np.array(vecs).T
    print >> out
    print >> out
    print >> out, 'augmented distance matrix:'
    print >> out, D_aug
    print >> out
    print >> out, 'scaled Fiedler vector of augmented Laplacian limit:'
    print >> out, fiedler_aug
    print >> out
    print >> out, 'eigenvalues of pinv of augmented Laplacian limit:'
    print >> out, vals_aug
    print >> out
    print >> out, 'rows are eigenvectors of pinv of augmented Laplacian limit:'
    print >> out, np.array(vecs_aug)
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
