"""
Test some 3D plotting with mayavi2.

The idea is that the mlab module tracks a bunch of state
as the script is executed,
and then at the end, the plot is shown in a little window.
"""

from enthought.mayavi import mlab
import numpy as np

from StringIO import StringIO
import math
from itertools import product

import NewickIO
import Euclid
import MatrixUtil
import EigUtil
import FelTree
import ProofDecoration
import const

g_tree_string = const.read('20100730g').rstrip()

g_plane_opacity = 0.1
g_crossing_opacity = 1.0
g_crossing_radius = 0.1


def add_yz_plane():
    x = [0, 0, 0, 0]
    y = [-1, -1, 1, 1]
    z = [-1, 1, 1, -1]
    tris = [(0, 3, 2), (2, 1, 0)]
    mesh = mlab.triangular_mesh(x, y, z,
            tris,
            color=(1, 0, 0),
            opacity=g_plane_opacity)

def add_zx_plane():
    y = [0, 0, 0, 0]
    x = [-1, -1, 1, 1]
    z = [-1, 1, 1, -1]
    tris = [(0, 3, 2), (2, 1, 0)]
    mesh = mlab.triangular_mesh(x, y, z,
            tris,
            color=(0, 1, 0),
            opacity=g_plane_opacity)

def add_xy_plane():
    z = [0, 0, 0, 0]
    x = [-1, -1, 1, 1]
    y = [-1, 1, 1, -1]
    tris = [(0, 3, 2), (2, 1, 0)]
    mesh = mlab.triangular_mesh(x, y, z,
            tris,
            color=(0, 0, 1),
            opacity=g_plane_opacity)

def draw_3d_tree(X, Y, Z, index_edges):
    for a, b in index_edges:
        foo = mlab.plot3d(
                [X[a], X[b]],
                [Y[a], Y[b]],
                [Z[a], Z[b]],
                color=(.5, .5, .5))

def gen_angles(n):
    for i in range(n+1):
        yield 2*math.pi*i/float(n)

def draw_x_crossing(y, z):
    angles = list(gen_angles(20))
    ys = [y + g_crossing_radius*math.cos(a) for a in angles]
    zs = [z + g_crossing_radius*math.sin(a) for a in angles]
    xs = np.zeros_like(ys)
    foo = mlab.plot3d(
            xs, ys, zs,
            color=(1, 0, 0),
            opacity=g_crossing_opacity)

def draw_y_crossing(x, z):
    angles = list(gen_angles(20))
    xs = [x + g_crossing_radius*math.cos(a) for a in angles]
    zs = [z + g_crossing_radius*math.sin(a) for a in angles]
    ys = np.zeros_like(xs)
    foo = mlab.plot3d(
            xs, ys, zs,
            color=(0, 1, 0),
            opacity=g_crossing_opacity)

def draw_z_crossing(x, y):
    angles = list(gen_angles(20))
    xs = [x + g_crossing_radius*math.cos(a) for a in angles]
    ys = [y + g_crossing_radius*math.sin(a) for a in angles]
    zs = np.zeros_like(xs)
    foo = mlab.plot3d(
            xs, ys, zs,
            color=(0, 0, 1),
            opacity=g_crossing_opacity)

def draw_crossings(X, Y, Z, index_edges):
    for a, b in index_edges:
        if X[a]*X[b] < 0:
            t = abs(X[b]) / abs(X[b] - X[a])
            draw_x_crossing(
                    t*Y[a] + (1-t)*Y[b],
                    t*Z[a] + (1-t)*Z[b])
        if Y[a]*Y[b] < 0:
            t = abs(Y[b]) / abs(Y[b] - Y[a])
            draw_y_crossing(
                    t*X[a] + (1-t)*X[b],
                    t*Z[a] + (1-t)*Z[b])
        if Z[a]*Z[b] < 0:
            t = abs(Z[b]) / abs(Z[b] - Z[a])
            draw_z_crossing(
                    t*X[a] + (1-t)*X[b],
                    t*Y[a] + (1-t)*Y[b])


def get_grant_proposal_points_b(lfdi):
    M, p, q = lfdi.M, lfdi.p, lfdi.q
    G = -.5 * M
    GQ, GX, GXT, GP = ProofDecoration.get_corners(G, q, p)
    # Get the eigendecomposition of the leaf-only Gower matrix.
    ws, vs = EigUtil.eigh(GQ)
    S = np.diag(ws)
    U = np.vstack(vs).T
    USUT = np.dot(np.dot(U, S), U.T)
    if not np.allclose(USUT, GQ):
        raise ValueError('eigenfail')
    S_sqrt = np.diag(np.sqrt(ws))
    X = np.dot(U, S_sqrt)
    # Find the imputed internal points.
    S_sqrt_pinv = np.linalg.pinv(S_sqrt)
    #W = np.dot(np.dot(S_sqrt_pinv, GX.T), U)
    try:
        W = np.dot(np.dot(GX.T, U), S_sqrt_pinv)
    except ValueError as e:
        arr = [
                GX.shape,
                U.shape,
                S_sqrt_pinv.shape]
        msg = ', '.join(str(x) for x in arr)
        raise ValueError(msg)
    # put them together and get only the first coordinates
    full_points = np.vstack([X, W])
    X = full_points.T[0]
    Y = full_points.T[1]
    Z = full_points.T[2]
    return X, Y, Z

def get_index_edges(tree, ordered_ids):
    """
    Given a tree and some ordered ids, get edges defined on indices.
    @param tree: the tree object
    @param ordered_ids: the returned index pairs are for this sequence
    @return: a collection of index pairs defining edges
    """
    # map ids to indices
    id_to_index = dict((myid, index) for index, myid in enumerate(ordered_ids))
    # each edge in this set is a frozenset of two indices
    index_edges = set()
    for node in tree.preorder():
        index = id_to_index[id(node)]
        for neighbor in node.gen_neighbors():
            neighbor_index = id_to_index[id(neighbor)] 
            index_edges.add(frozenset([index, neighbor_index]))
    return index_edges

def add_tree():
    # construct the matrices to be used for the eigendecomposition
    lfdo = ProofDecoration.tree_string_to_LFDO(g_tree_string)
    lfdi = ProofDecoration.LFDO_to_LFDI(lfdo)
    # we need the ordered ids themselves to to construct the edges
    tree = NewickIO.parse(g_tree_string, FelTree.NewickTree)
    ordered_ids = ProofDecoration.tree_to_leaf_first_ids(tree)
    index_edges = get_index_edges(tree, ordered_ids)
    # define the points
    X, Y, Z = get_grant_proposal_points_b(lfdi)
    # draw the image
    draw_3d_tree(X, Y, Z, index_edges)
    draw_crossings(X, Y, Z, index_edges)


# try to not pop up an annoying window
mlab.options.offscreen = True

add_yz_plane()
add_zx_plane()
add_xy_plane()
add_tree()

#mlab.show()
mlab.savefig('myfig.png', size=(640, 480))
#mlab.savefig('mymodel.vrml', size=(640, 480))
