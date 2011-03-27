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


g_line_width_thin = 1.0
g_line_width_normal = 2.0
g_line_width_thick = 4.0

g_color_background = (0.9, 0.9, 0.9)
g_color_dark = (0.0, 0.0, 0.0)
g_color_light = (0.6, 0.6, 0.6)
g_color_bad_edge = (1.0, 0, 0)

g_barb_radius = 5.0

g_min_valuation_radius = 1e-8


def wat(tree, eig_idx1, eig_idx2):
    """
    @param tree: a SpatialTree
    @param eig_idx1: first eigen index, 1 is Fiedler
    @param eig_idx2: second eigen index, 1 is Fiedler
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
    layout = FastDaylightLayout.StraightBranchLayout()
    layout.do_layout(tree)
    # draw the tree
    return get_tree_image(tree, (640, 480), ext, id_to_v1, id_to_v2)

def is_bad_edge(node, child, v1, v2):
    v1_pair = (v1[id(node)], v1[id(child)])
    v2_pair = (v2[id(node)], v2[id(child)])
    if max(abs(v) for v in v1_pair) < g_min_valuation_radius:
        return True
    if max(abs(v) for v in v2_pair) < g_min_valuation_radius:
        return True
    return False

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

