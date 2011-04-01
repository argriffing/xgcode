"""
Draw a tree annotated with roots of eigenfunctions.
"""

import math
import string

import cairo

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

g_border_outer = 15
g_border_inner = 40
g_min_pane_width = 10
g_min_pane_height = 10

g_vertex_dot_margin = 0.8
g_vertex_dot_radius_thin = 1.2
g_vertex_dot_radius_thick = 2

g_label_distance = 8


def is_bad_edge(node, child, v1, v2):
    v1_pair = (v1[id(node)], v1[id(child)])
    v2_pair = (v2[id(node)], v2[id(child)])
    if max(abs(v) for v in v1_pair) < g_min_valuation_radius:
        return True
    if max(abs(v) for v in v2_pair) < g_min_valuation_radius:
        return True
    return False

def get_eg2_image(tree, max_size, image_format, v1, v2s,
        bdrawbackground, bdrawvertices, bdrawlabels):
    """
    Get the image of the tree.
    This could be called from outside the module.
    @param tree: something like a SpatialTree
    @param max_size: (max_width, max_height)
    @param image_format: a string that determines the image format
    @param v1: a valuation dictionary
    @param v2s: fourteen valuation dictionaries
    @param bdrawbackground: flag to draw the background
    @param bdrawvertices: flag to draw vertices
    @param bdrawlabels: flag to draw labels
    @return: a string containing the image data
    """
    # Define the number of columns and the number of rows.
    # Note that there is raggedness.
    ncols = 3
    nrows = 5
    # hardcoded
    # TODO change the npairs name as it is not applicable
    npairs = 14
    # get the max width and height per pane
    max_width, max_height = max_size
    max_pane_width = (
            max_width - 2*g_border_outer - g_border_inner*(ncols-1))/ncols
    max_pane_height = (
            max_height - 2*g_border_outer - g_border_inner*(nrows-1))/nrows
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
    width = 2*g_border_outer + (ncols-1)*g_border_inner + ncols*pane_width
    height = 2*g_border_outer + (nrows-1)*g_border_inner + nrows*pane_height
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
        # draw rectangles for the panes
        if False:
            gap = g_border_inner - 2*g_border_outer
            if gap > 0:
                for row in range(nrows):
                    for col in range(ncols):
                        index = col*nrows + row
                        if index < npairs:
                            context.save()
                            context.set_source_rgb(*g_color_background)
                            x = col*(pane_width + g_border_inner)
                            y = row*(pane_height + g_border_inner)
                            w = 2*g_border_outer + pane_width
                            h = 2*g_border_outer + pane_height
                            context.rectangle(x, y, x+w, y+h)
                            context.fill()
                            context.restore()
    # draw the trees
    context.translate(g_border_outer + pane_width/2.0, 0)
    context.translate(0, g_border_outer + pane_height/2.0)
    for i, v2 in enumerate(v2s):
        # draw the tree into the context
        if bdrawbackground:
            bgcolor = g_color_background
        else:
            # Pretend that if there is not a background color
            # then the background color is white.
            # FIXME
            bgcolor = (1.0, 1.0, 1.0)
        draw_single_tree(tree, context, v1, v2, bgcolor,
                bdrawvertices, bdrawlabels)
        # draw the pane label into the context
        # TODO conditional
        if False:
            if i < len(string.uppercase):
                letter = string.uppercase[i]
            else:
                letter = '?'
            context.save()
            context.set_font_size(20.0)
            xbear, ybear, w, h, xadv, yadv = context.text_extents(letter)
            xtarget = -pane_width/2
            ytarget = -pane_height/2 + h
            context.move_to(xtarget, ytarget)
            context.show_text(letter)
            context.restore()
        # move the drawing position
        if i % nrows == nrows - 1:
            # move to the next column
            context.translate(g_border_inner + pane_width, 0)
            # move to the first row
            context.translate(0, -(nrows-1)*(g_border_inner + pane_height))
        else:
            context.translate(0, g_border_inner + pane_height)
    # get the image string
    return cairo_helper.get_image_string()

def get_forest_image(tree, max_size, image_format, vs, ncols,
        bdrawbackground, bdrawvertices, bdrawlabels):
    """
    Get the image of the tree.
    This could be called from outside the module.
    @param tree: something like a SpatialTree
    @param max_size: (max_width, max_height)
    @param image_format: a string that determines the image format
    @param vs: sequence of maps from node id to valuation
    @param ncols: use this many columns
    @param bdrawbackground: flag to draw the background
    @param bdrawvertices: flag to draw vertices
    @param bdrawlabels: flag to draw labels
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
    max_pane_width = (
            max_width - 2*g_border_outer - g_border_inner*(ncols-1))/ncols
    max_pane_height = (
            max_height - 2*g_border_outer - g_border_inner*(nrows-1))/nrows
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
    width = 2*g_border_outer + (ncols-1)*g_border_inner + ncols*pane_width
    height = 2*g_border_outer + (nrows-1)*g_border_inner + nrows*pane_height
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
        # draw rectangles for the panes
        if False:
            gap = g_border_inner - 2*g_border_outer
            if gap > 0:
                for row in range(nrows):
                    for col in range(ncols):
                        index = col*nrows + row
                        if index < npairs:
                            context.save()
                            context.set_source_rgb(*g_color_background)
                            x = col*(pane_width + g_border_inner)
                            y = row*(pane_height + g_border_inner)
                            w = 2*g_border_outer + pane_width
                            h = 2*g_border_outer + pane_height
                            context.rectangle(x, y, x+w, y+h)
                            context.fill()
                            context.restore()
    # draw the trees
    context.translate(g_border_outer + pane_width/2.0, 0)
    context.translate(0, g_border_outer + pane_height/2.0)
    for i, (v1, v2) in enumerate(iterutils.pairwise(vs)):
        # draw the tree into the context
        if bdrawbackground:
            bgcolor = g_color_background
        else:
            # Pretend that if there is not a background color
            # then the background color is white.
            # FIXME
            bgcolor = (1.0, 1.0, 1.0)
        draw_single_tree(tree, context, v1, v2, bgcolor,
                bdrawvertices, bdrawlabels)
        # draw the pane label into the context
        # TODO conditional
        if i < len(string.uppercase):
            letter = string.uppercase[i]
        else:
            letter = '?'
        context.save()
        context.set_font_size(20.0)
        xbear, ybear, w, h, xadv, yadv = context.text_extents(letter)
        xtarget = -pane_width/2
        ytarget = -pane_height/2 + h
        context.move_to(xtarget, ytarget)
        context.show_text(letter)
        context.restore()
        # move the drawing position
        if i % nrows == nrows - 1:
            # move to the next column
            context.translate(g_border_inner + pane_width, 0)
            # move to the first row
            context.translate(0, -(nrows-1)*(g_border_inner + pane_height))
        else:
            context.translate(0, g_border_inner + pane_height)
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
    bgcolor = g_color_background
    bdrawlabels = False
    bdrawvertices = False
    draw_single_tree(tree, context, v1, v2, bgcolor,
            bdrawlabels, bdrawvertices)
    # get the image string
    return cairo_helper.get_image_string()

def _draw_bad_branches(tree, context, v1, v2):
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

def _draw_directed_branches(tree, context, v1, v2):
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

def _draw_ticks(tree, context, v1, v2):
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

def get_free_angle(tree, node):
    """
    @param tree: something like a SpatialTree
    @param node: something like a SpatialNode
    @return: get the angle from the node to some free space
    """
    # FIXME this is implemented quadratic but it could be implemented linear.
    # get the list of angles away from the node of interest
    angles = []
    origin = tree._layout_to_display(node.location)
    for a, b in tree.gen_bidirected_branches():
        if a is node:
            target = tree._layout_to_display(b.location)
            theta = SpatialTree.get_angle(origin, target)
            angles.append(theta)
    # if there is only one angle then return its reflection
    if len(angles) == 1:
        return angles[0] + math.pi
    # If there are multiple angles then get the mid angle
    # of the widest interval.
    # Begin by sorting the angles.
    angles.sort()
    # convert pairs of sorted angles to angle intervals
    intervals = []
    for low, high in iterutils.pairwise(angles + [angles[0]]):
        intervals.append(SpatialTree.AngleInterval(low, high))
    # return the mid angle of the widest interval
    mag, interval = max((x.get_magnitude(), x) for x in intervals)
    return interval.get_mid_angle()

def draw_single_tree(tree, context, v1, v2, bgcolor,
        bdrawvertices, bdrawlabels):
    """
    This is most likely called only from inside the module.
    @param tree: a fitted SpatialTree
    @param context: a cairo context with origin at tree center
    @param v1: maps node id to valuation for line thickness
    @param v2: maps node id to valuation for zero crossing ticks
    @param bgcolor: background color
    @param bdrawvertices: flag to draw vertices
    @param bdrawlabels: flag to draw labels
    """
    _draw_bad_branches(tree, context, v1, v2)
    _draw_directed_branches(tree, context, v1, v2)
    if bdrawvertices:
        # make a set of nodes ids which are involved in bad edges
        bad_node_ids = set()
        for node, child in tree.gen_bidirected_branches():
            if is_bad_edge(node, child, v1, v2):
                bad_node_ids.add(id(node))
        # get (node, radius, color) triples
        triples = []
        for node in tree.preorder():
            d = id(node)
            if d in bad_node_ids:
                r = g_vertex_dot_radius_thick
                fgcolor = g_color_bad_edge
            elif v1[d] < 0:
                r = g_vertex_dot_radius_thin
                fgcolor = g_color_light
            else:
                r = g_vertex_dot_radius_thick
                fgcolor = g_color_light
            triples.append((node, r, fgcolor))
        # erase a disk around each vertex
        for node, radius, fgcolor in triples:
            # get the parameters for the blank dot
            r = radius + g_vertex_dot_margin
            x, y = tree._layout_to_display(node.location)
            # draw the blank dot
            context.save()
            context.set_source_rgb(*bgcolor)
            context.arc(x, y, r, 0, 2 * math.pi)
            context.fill()
            context.restore()
        # Draw a disk around each vertex
        # depending on its sign status and its bad edge status.
        for node, radius, fgcolor in triples:
            # get the parameters for the blank dot
            r = radius
            x, y = tree._layout_to_display(node.location)
            # draw the dot
            context.save()
            context.set_source_rgb(*fgcolor)
            context.arc(x, y, r, 0, 2 * math.pi)
            context.fill()
            context.restore()
    _draw_ticks(tree, context, v1, v2)
    if bdrawlabels:
        # draw the node labels
        for node in tree.preorder():
            label = node.get_name()
            if label:
                # get the parameters for the label
                theta = get_free_angle(tree, node)
                x, y = tree._layout_to_display(node.location)
                xlab = x + g_label_distance * math.cos(theta)
                ylab = y + g_label_distance * math.sin(theta)
                # draw the text, centered on the target point
                context.save()
                xbear, ybear, w, h, xadv, yadv = context.text_extents(label)
                xtarget = xlab - w/2
                ytarget = ylab + h/2
                context.move_to(xtarget, ytarget)
                context.show_text(label)
                context.restore()

