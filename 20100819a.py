"""Look at properties of the limit of the augmented distance matrix.
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
import RUtil

g_tree_string = '((1:1, 2:0.5)6:1, 5:1, (3:0.2, 4:0.5)7:1)8;'

def get_form():
    """
    @return: a list of form objects
    """
    # define the form objects
    form_objects = [
            Form.MultiLine('tree_string', 'newick tree',
                g_tree_string),
            Form.CheckGroup('options', 'output options', [
                Form.CheckItem('show_aug',
                    'show the augmented distance matrix')])]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_alpha_multiplier(D, ntips):
    n = len(D)
    q = ntips
    p = n - ntips
    Xt = D[p:].T[:p].T
    Dq = D[p:].T[p:].T
    uno = np.ones((q,p), dtype=float)
    right = np.dot(Dq, uno)/q
    R = Xt - right
    return R - np.mean(R, axis=0)

def get_beta_multiplier(D, ntips):
    n = len(D)
    q = ntips
    p = n - ntips
    X = D[p:].T[:p]
    Dq = D[p:].T[p:].T
    uno = np.ones((p,q), dtype=float)
    right = np.dot(uno, Dq)/q
    R = X - right
    return (R.T - np.mean(R.T, axis=0)).T

def get_response_content(fs):
    # build the newick tree from the string
    tree = NewickIO.parse(fs.tree_string, FelTree.NewickTree)
    nvertices = len(list(tree.preorder()))
    nleaves = len(list(tree.gen_tips()))
    ninternal = nvertices - nleaves
    # get ordered ids with the internal nodes first
    ordered_ids = get_ordered_ids(tree)
    leaf_ids = [id(node) for node in tree.gen_tips()]
    # get the distance matrix and the augmented distance matrix
    D_leaf = np.array(tree.get_partial_distance_matrix(leaf_ids))
    D = np.array(tree.get_partial_distance_matrix(ordered_ids))
    # analyze the leaf distance matrix
    X_leaf = Euclid.edm_to_points(D_leaf)
    w_leaf, v_leaf = EigUtil.eigh(Euclid.edm_to_dccov(D_leaf))
    V_leaf = np.array(v_leaf).T
    # explicitly compute the limiting points as the number of dups increases
    X = Euclid.edm_to_points(D)
    X -= np.mean(X[-nleaves:], axis=0)
    XL = X[-nleaves:]
    U, s, Vt = np.linalg.svd(XL)
    Z = np.dot(X, Vt.T)
    # hack the Z matrix to show the leaf-related eigenvectors
    Z = Z.T[:nleaves-1].T
    WY = Z / np.sqrt(w_leaf[:-1])
    # compute a product using the first few rows of WY
    W = WY[:ninternal]
    M_alpha = get_alpha_multiplier(D, nleaves)
    MW_alpha = np.dot(M_alpha, W)
    # compute a product using the first few rows of WY
    M_beta = get_beta_multiplier(D, nleaves)
    MY_beta = np.dot(M_beta, V_leaf)
    # report the results
    np.set_printoptions(linewidth=300, threshold=10000)
    out = StringIO()
    print >> out, 'leaf distance matrix:'
    print >> out, D_leaf
    print >> out
    print >> out, 'eigenvalues derived from the leaf distance matrix'
    print >> out, w_leaf
    print >> out
    print >> out, 'corresponding eigenvectors (as columns)'
    print >> out, V_leaf
    print >> out
    print >> out, "candidates for [W' Y']':"
    print >> out, WY
    print >> out
    print >> out, 'candidates for W:'
    print >> out, W
    print >> out
    print >> out, 'left multiplier of W:'
    print >> out, M_alpha
    print >> out
    print >> out, 'each column is a (left multiplier, W) product:'
    print >> out, MW_alpha
    print >> out
    print >> out, 'left multiplier of Y:'
    print >> out, M_beta
    print >> out
    print >> out, 'each column is a (left multiplier, Y) product:'
    print >> out, MY_beta
    print >> out
    print >> out, 'the above matrix divided by 2*eigenvalue:'
    print >> out, MY_beta / (2*np.array(w_leaf))
    print >> out
    return out.getvalue()

def get_ordered_ids(tree):
    """
    Maybe I could use postorder here instead.
    @param tree: a tree
    @return: a list of ids beginning with the leaves
    """
    ordered_ids = []
    ordered_ids.extend(id(node) for node in tree.gen_internal_nodes())
    ordered_ids.extend(id(node) for node in tree.gen_tips())
    return ordered_ids
