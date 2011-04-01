"""
Get a harmonic extension.
Use the harmonic extensions of eigenvectors
of the Schur complement Laplacian matrix.
"""

import numpy as np
import scipy


def get_harmonic_valuations(tree, eig_idx):
    """
    @param tree: a Newick tree
    @param eig_idx: eigen index, 1 is Fiedler
    @return: map from node id to harmonic valuation
    """
    # make the adjacency matrix
    ordered_tip_ids = [id(node) for node in tree.gen_tips()]
    ordered_internal_ids = [id(node) for node in tree.gen_internal_nodes()]
    ordered_ids = ordered_tip_ids + ordered_internal_ids
    id_to_idx = dict((myid, i) for i, myid in enumerate(ordered_ids))
    q = len(ordered_tip_ids)
    p = len(ordered_internal_ids)
    N = q + p
    A = np.zeros((N,N))
    for na, nb, blen in tree.gen_bidirected_branches_with_length():
        weight = 1/float(blen)
        idxa = id_to_idx[id(na)]
        idxb = id_to_idx[id(nb)]
        A[idxa, idxb] = weight
    # check the requested indices
    eig_msg = 'eigenfunction indices must be less than the number of leaves'
    if eig_idx >= q:
        raise ValueError(eig_msg)
    # define the Laplacian matrix and its pieces
    L = np.diag(np.sum(A, axis=0)) - A
    L11 = L[:q][:, :q]
    L12 = L[:q][:, -p:]
    L22 = L[-p:][:, -p:]
    L22_pinv = np.linalg.pinv(L22)
    L_star = L11 - np.dot(L12, np.dot(L22_pinv, L12.T))
    W, V1 = scipy.linalg.eigh(L_star)
    V2 = -np.dot(np.dot(L22_pinv, L12.T), V1)
    V = np.vstack([V1, V2])
    # define the vertex valuations
    id_to_v = dict((myid, V[i, eig_idx]) for i, myid in enumerate(
        ordered_ids))
    return id_to_v

