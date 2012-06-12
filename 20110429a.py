"""
Draw a combination of the 2D tip MDS and the harmonic extension of two trees.

The first input tree supplies the 2D MDS relating the leaves.
The second input tree supplies the topology and branch lengths
for the harmonic extension.
This gives a way to draw a proposed tree
inside a given 2D MDS.
This visualization could be generalized to the case where
the 2D MDS does not come from a tree,
but this generalization is not implemented here.
"""

from StringIO import StringIO

import numpy as np
import cairo
import scipy.linalg

import Form
import FormOut
import Ftree
import FtreeIO
import CairoUtil


def get_form():
    """
    @return: the body of a form
    """
    # define default tree strings
    true_s = '(b:2.622, d:1.013, (e:1.496, (a:2.749, c:0.338):1.474):0.889);'
    test_s = '(b:2.622, d:1.013, (e:1.496, (a:2.749, c:0.338):1.474):0.889);'
    # define the form objects
    form_objects = [
            Form.MultiLine('true_tree', 'leaf MDS tree', true_s),
            Form.MultiLine('test_tree', 'harmonic extension tree', test_s),
            Form.Float('scale', 'scale the image of the tree by this factor',
                200.0, low_exclusive=0.0),
            Form.ImageFormat()]
    return form_objects

def get_form_out():
    return FormOut.Image('mds')

def get_response_content(fs):
    # read the trees
    T_true, B_true, N_true = FtreeIO.newick_to_TBN(fs.true_tree)
    T_test, B_test, N_test = FtreeIO.newick_to_TBN(fs.test_tree)
    # we are concerned about the names of the leaves of the two trees
    true_leaves = Ftree.T_to_leaves(T_true)
    test_leaves = Ftree.T_to_leaves(T_test)
    true_leaf_to_n = dict((v, N_true[v]) for v in true_leaves)
    test_leaf_to_n = dict((v, N_test[v]) for v in test_leaves)
    # check that all leaves are named
    if len(true_leaves) != len(true_leaf_to_n):
        raise ValueError(
                'all leaves in the leaf MDS tree should be named')
    if len(test_leaves) != len(test_leaf_to_n):
        raise ValueError(
                'all leaves in the harmonic extension tree should be named')
    # check that within each tree all leaves are uniquely named
    if len(set(true_leaf_to_n.values())) != len(true_leaves):
        raise ValueError(
                'all leaf names in the leaf MDS tree should be unique')
    if len(set(test_leaf_to_n.values())) != len(test_leaves):
        raise ValueError(
                'all leaf names in the harmonic extension tree '
                'should be unique')
    # check that the leaf name sets are the same
    if set(true_leaf_to_n.values()) != set(test_leaf_to_n.values()):
        raise ValueError(
                'the two trees should have corresponding leaf names')
    # invert the leaf name maps
    true_n_to_leaf = dict((n, v) for v, n in true_leaf_to_n.items())
    test_n_to_leaf = dict((n, v) for v, n in test_leaf_to_n.items())
    # get correspondingly ordered leaf sequences
    leaf_names = true_leaf_to_n.values()
    true_leaves_reordered = [true_n_to_leaf[n] for n in leaf_names]
    test_leaves_reordered = [test_n_to_leaf[n] for n in leaf_names]
    # get the Schur complement matrix for the leaves
    L_schur_true = Ftree.TB_to_L_schur(T_true, B_true, true_leaves_reordered)
    # get the MDS points
    w, V = scipy.linalg.eigh(L_schur_true, eigvals=(1, 2))
    X = np.dot(V, np.diag(np.reciprocal(np.sqrt(w))))
    # get the linear operator that defines the harmonic extension
    test_internal = Ftree.T_to_internal_vertices(T_test)
    L22 = Ftree.TB_to_L_block(T_test, B_test,
            test_internal, test_internal)
    L21 = Ftree.TB_to_L_block(T_test, B_test,
            test_internal, test_leaves_reordered)
    M = -np.dot(np.linalg.pinv(L22), L21)
    # get the harmonic extension
    X_extension = np.dot(M, X)
    X_extended = np.vstack([X, X_extension])
    # draw the image
    v_to_index = Ftree.invseq(test_leaves_reordered + test_internal)
    physical_size = (640, 480)
    ext = Form.g_imageformat_to_ext[fs.imageformat]
    return get_animation_frame(ext, physical_size, fs.scale,
            v_to_index, T_test, X_extended)

def get_animation_frame(
        image_format, physical_size, scale,
        v_to_index, T, X):
    """
    This function is about drawing the tree.
    @param image_format: the image extension
    @param physical_size: the width and height of the image in pixels
    @param scale: a scaling factor
    @param v_to_index: maps vertices to their index
    @param T: defines the connectivity of the tree
    @param X: an array of 2D points
    @return: the animation frame as an image as a string
    """
    # before we begin drawing we need to create the cairo surface and context
    cairo_helper = CairoUtil.CairoHelper(image_format)
    surface = cairo_helper.create_surface(physical_size[0], physical_size[1])
    context = cairo.Context(surface)
    # define some helper variables
    x0 = physical_size[0] / 2.0
    y0 = physical_size[1] / 2.0
    npoints = len(X)
    # draw an off-white background
    context.save()
    context.set_source_rgb(.9, .9, .9)
    context.paint()
    context.restore()
    # draw the axes which are always in the center of the image
    context.save()
    context.set_source_rgb(.9, .7, .7)
    context.move_to(x0, 0)
    context.line_to(x0, physical_size[1])
    context.stroke()
    context.move_to(0, y0)
    context.line_to(physical_size[0], y0)
    context.stroke()
    context.restore()
    # draw the edges
    context.save()
    context.set_source_rgb(.8, .8, .8)
    for va, vb in T:
        a = v_to_index[va]
        b = v_to_index[vb]
        ax, ay = X[a].tolist()
        bx, by = X[b].tolist()
        context.move_to(x0 + ax*scale, y0 + ay*scale)
        context.line_to(x0 + bx*scale, y0 + by*scale)
        context.stroke()
    context.restore()
    # create the image
    return cairo_helper.get_image_string()
