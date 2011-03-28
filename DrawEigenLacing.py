"""
Draw a tree annotated with roots of eigenfunctions.
Use the harmonic extensions of eigenvectors
of the Schur complement Laplacian matrix.
"""

import math

import cairo
import numpy as np
import scipy

import Newick
import SpatialTree
import FastDaylightLayout
import CairoUtil
import iterutils


g_line_width_thin = 1.0
g_line_width_normal = 2.0
g_line_width_thick = 3.0

g_color_background = (0.9, 0.9, 0.9)
g_color_dark = (0.0, 0.0, 0.0)
g_color_light = (0.6, 0.6, 0.6)
g_color_bad_edge = (1.0, 0, 0)

g_barb_radius = 4.0

g_min_valuation_radius = 1e-8

g_pixel_border = 10
g_min_pane_width = 10
g_min_pane_height = 10


def get_harmonic_valuations(tree, eig_idx):
    """
    @param tree: a SpatialTree
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

def is_bad_edge(node, child, v1, v2):
    v1_pair = (v1[id(node)], v1[id(child)])
    v2_pair = (v2[id(node)], v2[id(child)])
    if max(abs(v) for v in v1_pair) < g_min_valuation_radius:
        return True
    if max(abs(v) for v in v2_pair) < g_min_valuation_radius:
        return True
    return False

def get_forest_image(tree, max_size, image_format, vs, ncols, bdrawbackground):
    """
    Get the image of the tree.
    This could be called from outside the module.
    @param tree: something like a SpatialTree
    @param max_size: (max_width, max_height)
    @param image_format: a string that determines the image format
    @param vs: sequence of maps from node id to valuation
    @param ncols: use this many columns
    @param bdrawbackground: whether or not to draw the background
    @return: a string containing the image data
    """
    npairs = len(vs)-1
    if npairs < 1:
        raise ValueError('not enough valuation maps')
    # get the number of rows
    nrows, remainder = divmod(npairs, ncols)
    if remainder:
        nrows += 1
    # get the max width and height per pane
    max_width, max_height = max_size
    max_pane_width = (max_width - g_pixel_border*(ncols+1))/ncols
    max_pane_height = (max_height - g_pixel_border*(nrows+1))/nrows
    if max_pane_width < g_min_pane_width:
        raise ValueError('not enough room')
    if max_pane_height < g_min_pane_height:
        raise ValueError('not enough room')
    max_pane_size = (max_pane_width, max_pane_height)
    # rotate and center the tree on (0, 0)
    tree.fit(max_pane_size)
    # get the width and height of the tree image
    xmin, ymin, xmax, ymax = tree.get_extents()
    pane_width = xmax - xmin
    pane_height = ymax - ymin
    width = g_pixel_border + ncols*(pane_width + g_pixel_border)
    height = g_pixel_border + nrows*(pane_height + g_pixel_border)
    # create the surface
    cairo_helper = CairoUtil.CairoHelper(image_format)
    surface = cairo_helper.create_surface(width, height)
    context = cairo.Context(surface)
    # draw the background
    if bdrawbackground:
        context.save()
        context.set_source_rgb(*g_color_background)
        context.paint()
        context.restore()
    # draw the trees
    context.translate(g_pixel_border + pane_width/2.0, 0)
    context.translate(0, g_pixel_border + pane_height/2.0)
    for i, (v1, v2) in enumerate(iterutils.pairwise(vs)):
        # draw the tree into the context
        draw_single_tree(tree, context, v1, v2)
        # move the drawing position
        if i % nrows == nrows - 1:
            # move to the next column
            context.translate(g_pixel_border + pane_width, 0)
            # move to the first row
            context.translate(0, -(nrows-1)*(g_pixel_border + pane_height))
        else:
            context.translate(0, g_pixel_border + pane_height)
    # get the image string
    return cairo_helper.get_image_string()

def get_single_tree_image(tree, max_size, image_format, v1, v2):
    """
    Get the image of the tree.
    This could be called from outside the module.
    @param tree: something like a SpatialTree
    @param max_size: (max_width, max_height)
    @param image_format: a string that determines the image format
    @param v1: maps node id to valuation for line thickness
    @param v2: maps node id to valuation for zero crossing ticks
    @return: a string containing the image data
    """
    # rotate and center the tree on (0, 0)
    tree.fit(max_size)
    # get the width and height of the tree image
    xmin, ymin, xmax, ymax = tree.get_extents()
    width = xmax - xmin
    height = ymax - ymin
    # create the surface
    cairo_helper = CairoUtil.CairoHelper(image_format)
    surface = cairo_helper.create_surface(width, height)
    context = cairo.Context(surface)
    # draw the background
    context.save()
    context.set_source_rgb(*g_color_background)
    context.paint()
    context.restore()
    # center on the tree
    context.translate(width/2.0, height/2.0)
    # draw the tree into the context
    draw_single_tree(tree, context, v1, v2)
    # get the image string
    return cairo_helper.get_image_string()

def draw_single_tree(tree, context, v1, v2):
    """
    This is most likely called only from inside the module.
    @param tree: a fitted SpatialTree
    @param context: a cairo context with origin at tree center
    @param v1: maps node id to valuation for line thickness
    @param v2: maps node id to valuation for zero crossing ticks
    """
    # draw the bad branches
    for node, child in tree.gen_directed_branches():
        if is_bad_edge(node, child, v1, v2):
            context.save()
            psrc = tree._layout_to_display(node.location)
            pdst = tree._layout_to_display(child.location)
            context.set_source_rgb(*g_color_bad_edge)
            context.move_to(*psrc)
            context.line_to(*pdst)
            context.stroke()
            context.restore()
    # draw the directed branches
    for node, child in tree.gen_directed_branches():
        if is_bad_edge(node, child, v1, v2):
            continue
        context.save()
        context.set_source_rgb(*g_color_light)
        # Get the valuations and (x,y) points of each node.
        vsrc, vdst = v1[id(node)], v1[id(child)]
        psrc = tree._layout_to_display(node.location)
        pdst = tree._layout_to_display(child.location)
        if vsrc < 0 and vdst < 0:
            context.set_line_width(g_line_width_thin)
            context.move_to(*psrc)
            context.line_to(*pdst)
            context.stroke()
        elif vsrc > 0 and vdst > 0:
            context.set_line_width(g_line_width_thick)
            context.move_to(*psrc)
            context.line_to(*pdst)
            context.stroke()
        else:
            # find the crossing point
            t = -vsrc / (vdst - vsrc)
            pmid = [psrc[i] + t*(pdst[i] - psrc[i]) for i in range(2)]
            # set line thickness for source
            if vsrc < 0:
                context.set_line_width(g_line_width_thin)
            else:
                context.set_line_width(g_line_width_thick)
            context.move_to(*psrc)
            context.line_to(*pmid)
            context.stroke()
            # set line thickness for destination
            if vdst < 0:
                context.set_line_width(g_line_width_thin)
            else:
                context.set_line_width(g_line_width_thick)
            context.move_to(*pmid)
            context.line_to(*pdst)
            context.stroke()
        context.restore()
    # draw the v2 zero crossing ticks perpendicular to the edge
    for node, child in tree.gen_directed_branches():
        if is_bad_edge(node, child, v1, v2):
            continue
        context.save()
        # Get the valuations and (x,y) points of each node.
        vsrc, vdst = v2[id(node)], v2[id(child)]
        psrc = tree._layout_to_display(node.location)
        pdst = tree._layout_to_display(child.location)
        if vsrc * vdst < 0:
            # find the crossing point
            t = -vsrc / (vdst - vsrc)
            pmid = [psrc[i] + t*(pdst[i] - psrc[i]) for i in range(2)]
            theta = math.atan2(pdst[1]-psrc[1], pdst[0]-psrc[0])
            barbx1 = pmid[0] + g_barb_radius * math.cos(theta + math.pi/2)
            barby1 = pmid[1] + g_barb_radius * math.sin(theta + math.pi/2)
            barbx2 = pmid[0] + g_barb_radius * math.cos(theta - math.pi/2)
            barby2 = pmid[1] + g_barb_radius * math.sin(theta - math.pi/2)
            # set line thickness for barb
            context.set_line_width(g_line_width_thin)
            context.move_to(barbx1, barby1)
            context.line_to(barbx2, barby2)
            context.stroke()
        context.restore()

