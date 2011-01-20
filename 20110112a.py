"""Draw an MDS with imputed internal nodes.

The purpose is to double check some equations in a manuscript.
"""


from StringIO import StringIO
import math

import numpy as np
import cairo

import Form
import FormOut
import NewickIO
import MatrixUtil
import EigUtil
import FelTree
import CairoUtil
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
    # construct the matrices to be used for the eigendecomposition
    lfdo = ProofDecoration.tree_string_to_LFDO(fs.tree_string)
    lfdi = ProofDecoration.LFDO_to_LFDI(lfdo)
    # we need the ordered ids themselves to to construct the edges
    tree = NewickIO.parse(fs.tree_string, FelTree.NewickTree)
    ordered_ids = ProofDecoration.tree_to_leaf_first_ids(tree)
    index_edges = get_index_edges(tree, ordered_ids)
    # define the points
    points = get_grant_proposal_points_b(lfdi)
    # draw the image
    ext = Form.g_imageformat_to_ext[fs.imageformat]
    return get_animation_frame(ext, physical_size, fs.scale,
            index_edges, points)

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
