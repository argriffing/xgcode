"""Draw an MDS with imputed internal nodes with the equations in the proposal.

This is a sanity check of the equations in the proposal.
"""


from StringIO import StringIO
import random
import os
import math
from itertools import product

import numpy as np
import cairo
import argparse

from SnippetUtil import HandlingError
import SnippetUtil
import Form
import FormOut
import NewickIO
import Euclid
import MatrixUtil
import EigUtil
import FelTree
import CairoUtil
import Progress
import ProofDecoration
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
            Form.Float('scale', 'scale the image of the tree by this factor',
                200.0, low_exclusive=0.0),
            Form.ImageFormat(),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.Image('tree')

def get_response_content(fs):
    # define the requested physical size of the images (in pixels)
    physical_size = (640, 480)
    # build the newick tree from the string
    tree = NewickIO.parse(fs.tree_string, FelTree.NewickTree)
    nvertices = len(list(tree.preorder()))
    nleaves = len(list(tree.gen_tips()))
    # Get ordered ids with the leaves first,
    # and get the corresponding distance matrix.
    lfdo = ProofDecoration.tree_string_to_LFDO(fs.tree_string)
    lfdi = ProofDecoration.LFDO_to_LFDI(lfdo)

    ordered_ids = get_ordered_ids_internal_first(tree)
    D = np.array(tree.get_partial_distance_matrix(ordered_ids))
    index_edges = get_index_edges(tree, ordered_ids)
    # draw the image
    ext = Form.g_imageformat_to_ext[fs.imageformat]
    points = get_grant_proposal_points(D, nleaves)
    return get_animation_frame(ext, physical_size, fs.scale,
            index_edges, points)

def get_grant_proposal_points(D, nleaves):
    """
    @return: rows are 2D points
    """
    B = D_to_B_third(D, nleaves)
    DQ = D_to_DQ_internal_first(D, nleaves)
    GQ = (-0.5) * MatrixUtil.double_centered(DQ)
    ws, vs = EigUtil.eigh(GQ)
    S = np.diag(ws)
    U = np.vstack(vs).T
    USUT = np.dot(np.dot(U, S), U.T)
    if not np.allclose(USUT, GQ):
        raise ValueError('eigenfail')
    S_sqrt = np.diag(np.sqrt(ws))
    S_sqrt_pinv = np.linalg.pinv(S_sqrt)
    # leaf-mds points
    X = np.dot(U, S_sqrt)
    # imputed internal points
    # the following line is as written in the proposal but it is wrong
    #W = np.dot(np.dot(S_sqrt_pinv, B), U)
    try:
        W = np.dot(np.dot(B, U), S_sqrt_pinv)
    except ValueError as e:
        arr = [
                B.shape,
                U.shape,
                S_sqrt_pinv.shape]
        msg = ', '.join(str(x) for x in arr)
        raise ValueError(msg)
    # put them together and get only the first coordinates
    full_points = np.vstack([W, X])
    points = full_points.T[:2].T
    return points

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

def get_animation_frame(
        image_format, physical_size, scale, index_edges, points):
    """
    This function is about drawing the tree.
    @param image_format: the image extension
    @param physical_size: the width and height of the image in pixels
    @param scale: a scaling factor
    @param index_edges: defines the connectivity of the tree
    @param points: an array of 2D points, the first few of which are leaves
    @return: the animation frame as an image as a string
    """
    # before we begin drawing we need to create the cairo surface and context
    cairo_helper = CairoUtil.CairoHelper(image_format)
    surface = cairo_helper.create_surface(physical_size[0], physical_size[1])
    context = cairo.Context(surface)
    # define some helper variables
    x0 = physical_size[0] / 2.0
    y0 = physical_size[1] / 2.0
    npoints = len(points)
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
    for edge in index_edges:
        ai, bi = tuple(edge)
        ax, ay = points[ai].tolist()
        bx, by = points[bi].tolist()
        context.move_to(x0 + ax*scale, y0 + ay*scale)
        context.line_to(x0 + bx*scale, y0 + by*scale)
        context.stroke()
    context.restore()
    # Draw vertices as translucent circles.
    context.save()
    context.set_source_rgba(0.2, 0.2, 1.0, 0.5)
    for point in points:
        x, y = point.tolist()
        nx = x0 + x*scale
        ny = y0 + y*scale
        dot_radius = 2.0
        context.arc(nx, ny, dot_radius, 0, 2*math.pi)
        context.fill()
    context.restore()
    # create the image
    return cairo_helper.get_image_string()

def get_ordered_ids_internal_first(tree):
    """
    @param tree: a tree
    @return: a list of ids
    """
    ordered_ids = []
    ordered_ids.extend(id(node) for node in tree.gen_internal_nodes())
    ordered_ids.extend(id(node) for node in tree.gen_tips())
    return ordered_ids

def D_to_DX_internal_first(D, q):
    """
    In the distance matrix the internal vertices come before the leaves
    @param D: the distance matrix with ordered rows and columns
    @param q: the number of leaves
    """
    p = len(D) - q
    return D[:p, p:]

def D_to_DQ_internal_first(D, q):
    """
    In the distance matrix the leaves come before the internal vertices
    @param D: the distance matrix with ordered rows and columns
    @param q: the number of leaves
    """
    p = len(D) - q
    return D[p:, p:]

def D_to_DX_leaves_first(D, q):
    """
    In the distance matrix the leaves come before the internal vertices
    @param D: the distance matrix with ordered rows and columns
    @param q: the number of leaves
    """
    return D[:q, q:]

def D_to_DU_leaves_first(D, q):
    """
    In the distance matrix the leaves come before the internal vertices
    @param D: the distance matrix with ordered rows and columns
    @param q: the number of leaves
    """
    return D[q:, q:]

def D_to_D_star(D, q):
    return D[:q, :q]

def get_column_mean_matrix(M):
    nrows = len(M)
    row = np.mean(M, axis=0)
    R = np.vstack([row]*M)
    if R.shape != M.shape:
        raise ValueError('internal shape fail')
    return R

def get_row_mean_matrix(M):
    R = get_column_mean_matrix(M.T).T
    if R.shape != M.shape:
        raise ValueError('internal shape fail')
    return R

def D_to_B_third(D, q):
    """
    This is the third attempt.
    It assumes that internal nodes are first.
    """
    p = len(D) - q
    DX = D_to_DX_internal_first(D, q)
    DQ = D_to_DQ_internal_first(D, q)
    DX_row_mean_matrix = get_row_mean_matrix(DX)
    DQ_column_mean_matrix = np.vstack([get_column_mean_matrix(DQ)[0]]*p)
    DQ_grand_mean = np.mean(DQ) * np.ones_like(DX)
    try:
        B = DX - DX_row_mean_matrix - DQ_column_mean_matrix + DQ_grand_mean
    except ValueError, e:
        arr = [
                DX.shape,
                DX_row_mean_matrix.shape,
                DQ_column_mean_matrix.shape,
                DQ_grand_mean.shape]
        msg = ', '.join(str(x) for x in arr)
        raise ValueError(msg)
    return B


def D_to_B(D, q):
    """
    This was the first attempt.
    It assumed that leaves were first.
    """
    DX = D_to_DX(D, q)
    DU = D_to_DU(D, q)
    DX_row_mean_matrix = get_row_mean_matrix(DX)
    DU_column_mean_matrix = np.vstack([get_column_mean_matrix(DU)[0]]*q)
    DU_grand_mean = np.mean(DU) * np.ones_like(DX)
    try:
        B = DX - DX_row_mean_matrix - DU_column_mean_matrix - DU_grand_mean
    except ValueError, e:
        arr = [
                DX.shape,
                DX_row_mean_matrix.shape,
                DU_column_mean_matrix.shape,
                DU_grand_mean.shape]
        msg = ', '.join(str(x) for x in arr)
        raise ValueError(msg)
    return B

def D_to_B_second(D, q):
    """
    I think this assumed that leaves were first
    """
    p = len(D)-q
    DX = D_to_DX(D, q)
    DS = D_to_D_star(D, q)
    row_mean_DS = np.mean(DS, axis=1)
    row_mean_matrix = np.vstack([row_mean_DS]*p).T
    col_mean_matrix = get_column_mean_matrix(DX)
    try:
        B = DX - col_mean_matrix - row_mean_matrix + np.mean(DS)
    except ValueError, e:
        arr = [
                DX.shape,
                col_mean_matrix.shape,
                row_mean_matrix.shape]
        msg = ', '.join(str(x) for x in arr)
        raise ValueError(msg)
    return B

def main(args):
    raise NotImplementedError('main does not do anything')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scale', type=float, default=1.0,
            help='define the drawing scale') 
    parser.add_argument('--physical_width', type=int, default=480,
            help='width (pixels)') 
    parser.add_argument('--physical_height', type=int, default=360,
            help='height (pixels)') 
    parser.add_argument('--tree', default=g_tree_string,
            help='newick tree with branch lengths')
    parser.add_argument('--image_format', default='png',
            choices=('png', 'svg', 'ps', 'pdf'),
            help='image format')
    parser.add_argument('--nframes', type=int, default=100,
            help='number of animation frames (image files) to create') 
    parser.add_argument('--interpolation', default='sigmoid',
            choices=('sigmoid', 'linear'),
            help='weights change according to this function')
    parser.add_argument('output_directory',
            help='path to the output directory for .png frames')
    main(parser.parse_args())

