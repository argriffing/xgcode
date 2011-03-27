""" Draw a tree annotated with roots of eigenfunctions.

Use the harmonic extensions of eigenvectors
of the Schur complement Laplacian matrix.
"""

import math

import cairo
import numpy as np
import scipy

from SnippetUtil import HandlingError
import SnippetUtil
import Newick
import SpatialTree
import FastDaylightLayout
import Form
import FormOut
import CairoUtil
import DrawEigenLacing


def get_form():
    """
    @return: the body of a form
    """
    # define the tree string
    tree_string = Newick.daylight_example_tree
    tree = Newick.parse(tree_string, Newick.NewickTree)
    formatted_tree_string = Newick.get_narrow_newick_string(tree, 60)
    # define the form objects
    form_objects = [
            Form.MultiLine('tree', 'newick tree', formatted_tree_string),
            Form.Integer('eig_idx1',
                'first eigenfunction index (1 means Fiedler)', 1, low=0),
            Form.Integer('eig_idx2',
                'second eigenfunction index (1 means Fiedler)', 2, low=0),
            Form.ImageFormat(),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.Image('tree')

def get_response_content(fs):
    # get a properly formatted newick tree with branch lengths
    tree = Newick.parse(fs.tree, SpatialTree.SpatialTree)
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
    if fs.eig_idx1 >= q or fs.eig_idx2 >= q:
        raise ValueError(eig_msg)
    # define the Laplacian matrix and its pieces
    L = np.diag(np.sum(A, axis=0)) - A
    L11 = L[:q][:, :q]
    L12 = L[:q][:, -p:]
    L22 = L[-p:][:, -p:]
    L22_pinv = np.linalg.pinv(L22)
    print L11.shape, L12.shape, L22.shape
    L_star = L11 - np.dot(L12, np.dot(L22_pinv, L12.T))
    W, V1 = scipy.linalg.eigh(L_star)
    V2 = -np.dot(np.dot(L22_pinv, L12.T), V1)
    V = np.vstack([V1, V2])
    # define the vertex valuations
    id_to_v1 = dict((myid, V[i, fs.eig_idx1]) for i, myid in enumerate(
        ordered_ids))
    id_to_v2 = dict((myid, V[i, fs.eig_idx2]) for i, myid in enumerate(
        ordered_ids))
    # do the layout
    try:
        layout = FastDaylightLayout.StraightBranchLayout()
        layout.do_layout(tree)
    except RuntimeError, e:
        pass
    # draw the image
    try:
        ext = Form.g_imageformat_to_ext[fs.imageformat]
        return DrawEigenLacing.get_single_tree_image(
                tree, (640, 480), ext, id_to_v1, id_to_v2)
    except CairoUtil.CairoUtilError, e:
        raise HandlingError(e)
