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
import SchurAlgebra
import iterutils


g_line_width_thin = 1.0
g_line_width_normal = 2.0
g_line_width_thick = 4.0

g_color_dark = (0.0, 0.0, 0.0)
g_color_light = (0.4, 0.4, 0.4)

g_barb_radius = 5.0


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
                'first eigenfunction index (1 means Fiedler)', 1, low=1),
            Form.Integer('eig_idx2',
                'second eigenfunction index (1 means Fiedler)', 2, low=1),
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
        return get_tree_image(tree, (640, 480), ext, id_to_v1, id_to_v2)
    except CairoUtil.CairoUtilError, e:
        raise HandlingError(e)


def get_tree_image(tree, max_size, image_format, v1, v2):
    """
    Get the image of the tree.
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
    context.set_source_rgb(.9, .9, .9)
    context.paint()
    context.restore()
    # center on the tree
    context.translate(width/2.0, height/2.0)
    # draw the directed branches
    for node, child in tree.gen_directed_branches():
        context.save()
        context.set_source_rgb(*g_color_light)
        # Get the valuations and (x,y) points of each node.
        vsrc = v1[id(node)]
        vdst = v1[id(child)]
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
        context.save()
        # Get the valuations and (x,y) points of each node.
        vsrc = v2[id(node)]
        vdst = v2[id(child)]
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
    # get the image string
    return cairo_helper.get_image_string()

