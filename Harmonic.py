"""
Get a harmonic extension.
Use the harmonic extensions of eigenvectors
of the Schur complement Laplacian matrix.
Note that scipy symmetric eigendecomposition
is a bit different from numpy symmetric eigendecomposition.
There might be differences in the ordering of the eigenpairs
and in the arbitrary sign assignments to the eigenvectors.
"""

import numpy as np
import scipy


def get_harmonic_valuations(tree, eig_idx):
    """
    @param tree: a Newick tree
    @param eig_idx: eigen index, 1 is Fiedler
    @return: map from node id to harmonic valuation
    """
    w, id_to_v = get_eigenvalue_and_harmonic_valuations(tree, eig_idx)
    return id_to_v

def get_eigenvalue_and_harmonic_valuations(tree, eig_idx):
    """
    @param tree: a Newick tree
    @param eig_idx: eigen index, 1 is Fiedler
    @return: eigenvalue, map from node id to harmonic valuation
    """
    # make the adjacency matrix
    ordered_tip_ids = [id(node) for node in tree.gen_tips()]
    ordered_internal_ids = [id(node) for node in tree.gen_internal_nodes()]
    ordered_ids = ordered_tip_ids + ordered_internal_ids
    id_to_idx = dict((myid, i) for i, myid in enumerate(ordered_ids))
    q = len(ordered_tip_ids)
    p = len(ordered_internal_ids)
    N = q + p
    # check the requested indices
    eig_msg = 'eigenfunction indices must be less than the number of leaves'
    if eig_idx >= q:
        raise ValueError(eig_msg)
    # define the Laplacian matrix and its pieces
    L = get_laplacian(tree, id_to_idx, q, p)
    L11 = L[:q][:, :q]
    L12 = L[:q][:, -p:]
    L22 = L[-p:][:, -p:]
    L22_pinv = np.linalg.pinv(L22)
    # get the eigenvectors and harmonic extensions of the Schur complement
    L_star = L11 - np.dot(L12, np.dot(L22_pinv, L12.T))
    W, V1 = scipy.linalg.eigh(L_star)
    V2 = -np.dot(np.dot(L22_pinv, L12.T), V1)
    V = np.vstack([V1, V2])
    # define the vertex valuations
    id_to_v = dict((myid, V[i, eig_idx]) for i, myid in enumerate(
        ordered_ids))
    return W[eig_idx], id_to_v

def get_laplacian(tree, id_to_idx, q, p):
    """
    @param tree: a Newick tree
    """
    N = q + p
    A = np.zeros((N,N))
    for na, nb, blen in tree.gen_bidirected_branches_with_length():
        weight = 1/float(blen)
        idxa = id_to_idx[id(na)]
        idxb = id_to_idx[id(nb)]
        A[idxa, idxb] = weight
    L = np.diag(np.sum(A, axis=0)) - A
    return L

def leaf_ids_to_schur_complement(tree, leaf_ids):
    """
    @param tree: a Newick tree
    @param leaf_ids: a sequence of leaf ids
    @return: the Schur complement matrix
    """
    internal_ids = [id(node) for node in tree.gen_internal_nodes()]
    ordered_ids = list(leaf_ids) + list(internal_ids)
    q = len(leaf_ids)
    p = len(internal_ids)
    id_to_idx = dict((myid, i) for i, myid in enumerate(ordered_ids))
    # define the Laplacian matrix and its pieces
    L = get_laplacian(tree, id_to_idx, q, p)
    L11 = L[:q][:, :q]
    L12 = L[:q][:, -p:]
    L22 = L[-p:][:, -p:]
    L22_pinv = np.linalg.pinv(L22)
    L_star = L11 - np.dot(L12, np.dot(L22_pinv, L12.T))
    return L_star

def leaf_names_to_schur_complement(tree, leaf_names):
    """
    @param tree: a Newick tree
    @param leaf_names: a sequence of leaf names
    @return: the Schur complement matrix
    """
    leaf_name_to_id = dict((node.name, id(node)) for node in tree.gen_tips())
    if set(leaf_names) != set(leaf_name_to_id):
        raise ValueError('leaf names do not match')
    ordered_leaf_ids = [leaf_name_to_id[name] for name in leaf_names]
    return leaf_ids_to_schur_complement(tree, leaf_ids)
