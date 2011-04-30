"""Draw a tree into the 2D MDS of an arbitrary leaf distance matrix.

The input distance matrix supplies the 2D MDS relating the leaves.
The input tree supplies the topology and branch lengths
for the harmonic extension.
This gives a way to draw a proposed tree
inside a given 2D MDS.
The default data is from
http://www.amjbot.org/cgi/content/full/93/9/1274
Genetic diversity and population structure in natural populations
of Moroccan Atlas cedar (Cedrus atlantica; Pinaceae)
determined with cpSSR markers.
The distance matrix is the Ds distance matrix
(not the Goldstein distance matrix).
The tree is from neighbor joining of the Goldstein distance matrix
with a few analysis peculiarities.
The first peculiarity is that a negative distance was forced to zero
before the neighbor joining analysis.
The second peculiarity is that negative branch lengths in the inferred
neighbor joining tree were transformed by taking their absolute values.
Furthermore the MDS (PCoA) as implemented by the NTSYS software
used by the authors of the reference paper ignores negative
eigenvalues of the distance matrix.
"Of course, an arbitrary dissimilarity matrix
may not be very compatible with a Euclidean metric.
In such cases many of the eigenvalues may be negative.
In performing such an analysis one hopes that such negative eigenvalues
are small and can be ignored."
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
import MatrixUtil
import Newick
import Util

# the default tree was created from this matrix
g_tree_generation_D_lines = [
        '0.00 0.21 1.19 0.79 0.18 1.57',
        '0.21 0.00 1.22 0.82 0.26 0.84',
        '1.19 1.22 0.00 -0.03 1.05 1.09',
        '0.79 0.82 -0.03 0.00 0.68 0.88',
        '0.18 0.26 1.05 0.68 0.00 0.64',
        '1.57 0.84 1.09 0.88 0.64 0.00']

g_default_D_lines = [
        '0.000 0.175 0.657 0.174 0.083 0.645',
        '0.175 0.000 0.572 0.200 0.144 0.477',
        '0.657 0.572 0.000 0.107 0.286 0.144',
        '0.174 0.200 0.107 0.000 0.119 0.314',
        '0.083 0.144 0.286 0.119 0.000 0.171',
        '0.645 0.477 0.144 0.314 0.171 0.000']

g_default_leaf_names = ['1', '2', '3', '4', '5', '6']

g_default_tree = '(5:0.06625, (6:0.521666666667, (3:0.1725, 4:0.1725):0.463333333333):0.32625, (1:0.1725, 2:0.0375):0.18125);'


def get_form():
    """
    @return: the body of a form
    """
    # format the default tree string
    tree = Newick.parse(g_default_tree, Newick.NewickTree)
    formatted_tree_string = Newick.get_narrow_newick_string(tree, 60)
    # make a real distance matrix
    D_s = [line.split() for line in g_default_D_lines]
    D_default = np.array([[float(x) for x in row] for row in D_s])
    # define the form objects
    form_objects = [
            Form.MultiLine('test_tree', 'harmonic extension tree',
                formatted_tree_string),
            Form.Matrix('D', 'leaf distance matrix',
                D_default, MatrixUtil.assert_symmetric),
            Form.MultiLine('names', 'ordered leaf names',
                '\n'.join(g_default_leaf_names)),
            Form.Float('scale', 'scale the image of the tree by this factor',
                200.0, low_exclusive=0.0),
            Form.ImageFormat(),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.Image('mds')

def MDS_v1(D):
    """
    Use the second and third eigenvalues regardless of negative eigenvalues.
    This method is probably not so great.
    @param D: tree distance matrix
    @return: leaf points
    """
    # get the Schur complement matrix analog for the leaves
    G = -0.5 * MatrixUtil.double_centered(D)
    L_schur_true = np.linalg.pinv(G)
    # get the MDS points
    w_all, V_all = scipy.linalg.eigh(L_schur_true)
    # get the relevant eigenvalues and eigenvectors
    imin = 0
    w = w_all[imin:imin+2]
    V = V_all[:, imin:imin+2]
    # get the MDS points
    return np.dot(V, np.diag(np.reciprocal(np.sqrt(w))))

def MDS_v2(D):
    """
    Use the two smallest positive eigenvalues.
    This method is the default for R PCoA software.
    @param D: tree distance matrix
    @return: leaf points
    """
    # get the Schur complement matrix analog for the leaves
    G = -0.5 * MatrixUtil.double_centered(D)
    L_schur_true = np.linalg.pinv(G)
    # get the MDS points
    w_all, V_all = scipy.linalg.eigh(L_schur_true)
    # get the index of the smallest positive eigenvalue
    eps = 1e-8
    xmin, imin = min((x, i) for i, x in enumerate(w_all) if x > eps)
    # get the relevant eigenvalues and eigenvectors
    w = w_all[imin:imin+2]
    V = V_all[:, imin:imin+2]
    # get the MDS points
    return np.dot(V, np.diag(np.reciprocal(np.sqrt(w))))

def MDS_v3(D):
    """
    Use the two smallest nonzero absolute value eigenvalues.
    This method is probably not so great.
    @param D: tree distance matrix
    @return: leaf points
    """
    # get the Schur complement matrix analog for the leaves
    G = -0.5 * MatrixUtil.double_centered(D)
    L_schur_true = np.linalg.pinv(G)
    # get the MDS points
    w_all, V_all = scipy.linalg.eigh(L_schur_true)
    # get the index of the smallest positive eigenvalue
    eps = 1e-8
    wi_sorted = sorted((abs(x), i) for i, x in enumerate(w_all))
    (xa, ia), (xb, ib) = wi_sorted[1:3]
    # get the relevant eigenvalues and eigenvectors
    w = w_all[[ia, ib]]
    V = V_all[:, [ia, ib]]
    # get the MDS points
    return np.dot(V, np.diag(np.reciprocal(np.sqrt(w))))

def MDS_v4(D):
    """
    Use the Lingoes method of removing negative eigenvalues.
    @param D: tree distance matrix
    @return: leaf points
    """
    eps = 1e-8
    # get the Schur complement matrix analog for the leaves
    G = -0.5 * MatrixUtil.double_centered(D)
    # do the Lingoes correction
    most_neg = scipy.linalg.eigh(G, eigvals_only=True, eigvals=(0,0))
    if most_neg < -eps:
        D_corrected = D - 2*(-most_neg)*(np.ones_like(D) - np.eye(len(D)))
        return MDS_v2(D_corrected)
    else:
        return MDS_v2(D)

def get_response_content(fs):
    # read the ordered leaf names for the distance matrix
    D_names = Util.get_stripped_lines(fs.names.splitlines())
    # read the tree
    T_test, B_test, N_test = FtreeIO.newick_to_TBN(fs.test_tree)
    # we are concerned about the names of the leaves of the two trees
    test_leaves = Ftree.T_to_leaves(T_test)
    test_leaf_to_n = dict((v, N_test[v]) for v in test_leaves)
    # check that all leaves are named
    if len(D_names) != len(fs.D):
        msg_a = 'the number of ordered leaf names should be the same '
        msg_b = 'as the number of rows in the distance matrix'
        raise HandlingError(msg_a + msg_b)
    if len(test_leaves) != len(test_leaf_to_n):
        msg = 'all leaves in the harmonic extension tree should be named'
        raise ValueError(msg)
    # check that leaves are uniquely named
    if len(set(D_names)) != len(D_names):
        msg = 'all ordered leaf names in the distance matrix should be unique'
        raise ValueError(msg)
    # check that the leaf name sets are the same
    if set(D_names) != set(test_leaf_to_n.values()):
        msg_a = 'the set of leaf names on the tree '
        msg_b = 'should be the same as '
        msg_c = 'the set of leaf names for the distance matrix'
        raise ValueError(msg_a + msg_b + msg_c)
    # invert the leaf name map
    test_n_to_leaf = dict((n, v) for v, n in test_leaf_to_n.items())
    # get correspondingly ordered leaf sequences
    test_leaves_reordered = [test_n_to_leaf[n] for n in D_names]
    # get the MDS points
    X = MDS_v4(fs.D)
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
