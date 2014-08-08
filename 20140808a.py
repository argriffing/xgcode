"""
Check vertex weight limit definition versus normalized Laplacian definition.

The normalized Laplacian definition uses the notation of Fan Chung.
The other vertex weight definition uses the notation of weighted MDS.
The pseudoinverse of the normalized Laplacian is related to a weighted
MDS matrix in a way that is similar to the way that the pseudoinverse
of a combinatorial Laplacian is related to the Gower matrix.

One question that I have not resolved is the calculation of the
appropriate MDS weights (or the Laplacian out-degrees) using the
original distance matrix and without going through an expensive
pseudoinversion.

"""
from StringIO import StringIO
import random
import math

import numpy as np
import scipy.linalg

import Form
import FormOut
import const
import NewickIO
import Euclid
import FelTree

g_tree_string = const.read('20100730g').rstrip()


def get_form():
    """
    @return: a list of form objects
    """
    # define the form objects
    form_objects = [
            Form.MultiLine('tree_string', 'newick tree', g_tree_string),
            ]
    return form_objects


def get_form_out():
    return FormOut.Report()


def get_ordered_ids(tree):
    """
    @param tree: a tree
    @return: a list of ids beginning with the leaves
    """
    ordered_ids = []
    ordered_ids.extend(id(node) for node in tree.gen_tips())
    ordered_ids.extend(id(node) for node in tree.gen_internal_nodes())
    return ordered_ids


def get_spectral_info(S):
    w = scipy.linalg.eigvalsh(S)
    eps = 1e-10
    wpinv = np.array([1/x if x > eps else 0 for x in w])
    return w, wpinv


def get_response_content(fs):

    # set up print options
    np.set_printoptions(
            linewidth=1000000,
            threshold=1000000,
            )
    out = StringIO()
    
    tree = NewickIO.parse(fs.tree_string, FelTree.NewickTree)

    # Get ordered ids with the leaves first.
    nvertices = len(list(tree.preorder()))
    nleaves = len(list(tree.gen_tips()))
    ordered_ids = get_ordered_ids(tree)

    # Report the full distance matrix.
    D_full = np.array(tree.get_partial_distance_matrix(ordered_ids))
    print >> out, 'full distance matrix:'
    print >> out, D_full
    print >> out

    # Extract the part of the distance matrix that relates only leaves.
    D = D_full[:nleaves, :nleaves]
    print >> out, 'leaf distance matrix:'
    print >> out, D
    print >> out

    # Compute the corresponding Laplacian matrix.
    L_comb = Euclid.edm_to_laplacian(D)
    w, wpinv = get_spectral_info(L_comb)
    print >> out, 'leaf combinatorial Laplacian matrix:'
    print >> out, L_comb
    print >> out, 'spectrum:', w
    print >> out, 'pinv spectrum:', wpinv
    print >> out

    # Compute the normalized Laplacian matrix.
    out_degrees = np.diag(L_comb)
    v = np.reciprocal(np.sqrt(out_degrees))
    L_norm = L_comb * np.outer(v, v)
    w, wpinv = get_spectral_info(L_norm)
    w, v = scipy.linalg.eigh(L_norm)
    print >> out, 'leaf normalized Laplacian matrix:'
    print >> out, L_norm
    print >> out, 'spectrum:', w
    print >> out, 'pinv spectrum:', wpinv
    print >> out, 'eigenvectors:'
    print >> out, v
    print >> out

    # Attempt to compute something related to weighted MDS.
    m = out_degrees
    M = np.diag(np.sqrt(m))
    I = np.identity(nleaves)
    e = np.ones(nleaves)
    E = I - np.outer(e, m) / np.inner(m, e)
    ME = np.dot(M, E)
    Q = -0.5 * ME.dot(D).dot(ME.T)
    w, wpinv = get_spectral_info(Q)
    w, v = scipy.linalg.eigh(Q)
    print >> out, 'a matrix related to weighted MDS:'
    print >> out, Q
    print >> out, 'spectrum:', w
    print >> out, 'pinv spectrum:', wpinv
    print >> out, 'eigenvectors:'
    print >> out, v
    print >> out


    # show the result
    return out.getvalue()

