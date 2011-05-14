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
import layout


g_line_width_thin = 1.0
g_line_width_normal = 2.0
g_line_width_thick = 3.0

g_color_background = (0.9, 0.9, 0.9)
g_color_dark = (0.0, 0.0, 0.0)
g_color_light = (0.6, 0.6, 0.6)
g_color_bad_edge = (1.0, 0, 0)
g_color_wavy_edge = (0.6, 0.6, 0.6)

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

g_wavelength = 3


def draw_wavy_line(ctx, x1, y1, x2, y2, wavelength):
    """
    @param ctx: cairo context
    @param x1: x coordinate of first point
    @param y1: y coordinate of first point
    @param x2: x coordinate of second point
    @param y2: y coordinate of second point
    @param wavelength: a value proportional to the wavelength
    """
    # set up the parameters
    d = math.hypot(y2-y1, x2-x1)
    nknots = 2 + int(math.ceil(d / wavelength))
    xincr = (x2-x1) / (nknots - 1)
    yincr = (y2-y1) / (nknots - 1)
    amplitude = math.hypot(xincr, yincr) / 2
    # start drawing
    ctx.save()
    ctx.set_line_width(g_line_width_thin)
    ctx.set_source_rgb(*g_color_wavy_edge)
    for i in range(nknots - 1):
        xa = x1 + i*xincr
        ya = y1 + i*yincr
        xb = xa + xincr
        yb = ya + yincr
        xm = (xa + xb)/2
        ym = (ya + yb)/2
        theta = math.atan2(yb-ya, xb-xa)
        if i % 2:
            theta += math.pi / 2
        else:
            theta -= math.pi / 2
        xc = xm + amplitude * math.cos(theta)
        yc = ym + amplitude * math.sin(theta)
        ctx.move_to(xa, ya)
        ctx.curve_to(xc, yc, xc, yc, xb, yb)
        ctx.stroke()
    ctx.restore()

def is_bad_edge_revised(node, child, v1):
    v1_pair = (v1[id(node)], v1[id(child)])
    if max(abs(v) for v in v1_pair) < g_min_valuation_radius:
        return True
    return False

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
        #FIXME
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
        # FIXME
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

def locations_to_extents(locations):
    """
    @param locations: (x, y) pairs
    @return: (xmin, ymin, xmax, ymax)
    """
    xmin = min(x for x, y in locations)
    ymin = min(y for x, y in locations)
    xmax = max(x for x, y in locations)
    ymax = max(y for x, y in locations)
    return (xmin, ymin, xmax, ymax)

def get_forest_image_revised(tree, max_size, image_format, vs,
        bdrawbackground, bdrawlabels, inner_margin, outer_margin,
        reflect_trees):
    """
    Get the image of the tree.
    This could be called from outside the module.
    The rectangle size is determined automatically
    so as to maximize the scaling factors of the trees in the panels.
    @param tree: something like a SpatialTree
    @param max_size: (max_width, max_height)
    @param image_format: a string that determines the image format
    @param vs: sequence of maps from node id to valuation
    @param bdrawbackground: flag to draw the background
    @param bdrawlabels: flag to draw labels
    @return: a string containing the image data
    """
    npairs = len(vs)
    if npairs < 1:
        raise ValueError('not enough valuation maps')
    # get all minimum rectangle sizes
    rect_sizes = layout.get_rect_sizes(npairs)
    # get (scaling_factor, w, h) triples
    triples = []
    max_width, max_height = max_size
    for ncols, nrows in rect_sizes:
        max_pane_width = (
                max_width - 2*outer_margin - inner_margin*(ncols-1))/ncols
        max_pane_height = (
                max_height - 2*outer_margin - inner_margin*(nrows-1))/nrows
        # require a minimum size
        if max_pane_width < g_min_pane_width:
            continue
        elif max_pane_height < g_min_pane_height:
            continue
        # append a triple after finding the scaling factor
        max_pane_size = (max_pane_width, max_pane_height)
        tree.fit(max_pane_size)
        triples.append((tree.scale, ncols, nrows))
    # if no triples were found then we fail
    if not triples:
        raise ValueError('not enough room')
    # get the nrows and ncols corresponding to the best scaling factor
    max_scale, ncols, nrows = max(triples)
    # get the position of each pane
    row_col_pairs = layout.min_rect_to_row_major(ncols, nrows, npairs)
    # re-fit the tree using the best size
    max_pane_width = (
            max_width - 2*outer_margin - inner_margin*(ncols-1))/ncols
    max_pane_height = (
            max_height - 2*outer_margin - inner_margin*(nrows-1))/nrows
    max_pane_size = (max_pane_width, max_pane_height)
    tree.fit(max_pane_size)
    # get the map from id to location for the final tree layout
    nodes = list(tree.preorder())
    id_to_location = dict(
            (id(n), tree._layout_to_display(n.location)) for n in nodes)
    if reflect_trees:
        id_to_location = dict(
                (v, (-x, y)) for v, (x,y) in id_to_location.items())
    # get the width and height of the tree image
    xmin, ymin, xmax, ymax = locations_to_extents(id_to_location.values())
    pane_width = xmax - xmin
    pane_height = ymax - ymin
    width = 2*outer_margin + (ncols-1)*inner_margin + ncols*pane_width
    height = 2*outer_margin + (nrows-1)*inner_margin + nrows*pane_height
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
    for i, ((row, col), (v1, v2)) in enumerate(
            zip(row_col_pairs, iterutils.pairwise(vs+[None]))):
        context.save()
        # center the context on the correct pane
        xtrans_initial = outer_margin + pane_width/2.0
        xtrans_extra = (inner_margin + pane_width)*col
        ytrans_initial = outer_margin + pane_height/2.0
        ytrans_extra = (inner_margin + pane_height)*row
        context.translate(
                xtrans_initial + xtrans_extra, ytrans_initial + ytrans_extra)
        # draw the tree into the context
        if bdrawbackground:
            bgcolor = g_color_background
        else:
            # Pretend that if there is not a background color
            # then the background color is white.
            bgcolor = (1.0, 1.0, 1.0)
        draw_single_tree_revised(
                tree, context, v1, v2, bgcolor, bdrawlabels, id_to_location)
        # draw the pane label into the context
        if i < len(string.uppercase):
            pane_label = str(i+1)
        else:
            pane_label = '?'
        context.save()
        #context.set_font_size(20.0)
        context.set_font_size(14.0)
        xbear, ybear, w, h, xadv, yadv = context.text_extents(pane_label)
        xtarget = -pane_width/2
        ytarget = -pane_height/2 + h
        context.move_to(xtarget, ytarget)
        context.show_text(pane_label)
        context.restore()
        # restore the position
        context.restore()
    # get the image string
    return cairo_helper.get_image_string()

def get_forest_image(tree, max_size, image_format, vs,
        bdrawbackground, bdrawvertices, bdrawlabels):
    """
    Get the image of the tree.
    This could be called from outside the module.
    The rectangle size is determined automatically
    so as to maximize the scaling factors of the trees in the panels.
    @param tree: something like a SpatialTree
    @param max_size: (max_width, max_height)
    @param image_format: a string that determines the image format
    @param vs: sequence of maps from node id to valuation
    @param bdrawbackground: flag to draw the background
    @param bdrawvertices: flag to draw vertices
    @param bdrawlabels: flag to draw labels
    @return: a string containing the image data
    """
    npairs = len(vs)-1
    if npairs < 1:
        raise ValueError('not enough valuation maps')
    # get all minimum rectangle sizes
    rect_sizes = layout.get_rect_sizes(npairs)
    # get (scaling_factor, w, h) triples
    triples = []
    max_width, max_height = max_size
    for ncols, nrows in rect_sizes:
        max_pane_width = (
                max_width - 2*g_border_outer - g_border_inner*(ncols-1))/ncols
        max_pane_height = (
                max_height - 2*g_border_outer - g_border_inner*(nrows-1))/nrows
        # require a minimum size
        if max_pane_width < g_min_pane_width:
            continue
        elif max_pane_height < g_min_pane_height:
            continue
        # append a triple after finding the scaling factor
        max_pane_size = (max_pane_width, max_pane_height)
        tree.fit(max_pane_size)
        triples.append((tree.scale, ncols, nrows))
    # if no triples were found then we fail
    if not triples:
        raise ValueError('not enough room')
    # get the nrows and ncols corresponding to the best scaling factor
    max_scale, ncols, nrows = max(triples)
    # get the position of each pane
    row_col_pairs = layout.min_rect_to_row_major(ncols, nrows, npairs)
    # re-fit the tree using the best size
    max_pane_width = (
            max_width - 2*g_border_outer - g_border_inner*(ncols-1))/ncols
    max_pane_height = (
            max_height - 2*g_border_outer - g_border_inner*(nrows-1))/nrows
    max_pane_size = (max_pane_width, max_pane_height)
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
    # draw the trees
    for i, ((row, col), (v1, v2)) in enumerate(
            zip(row_col_pairs, iterutils.pairwise(vs))):
        context.save()
        # center the context on the correct pane
        xtrans_initial = g_border_outer + pane_width/2.0
        xtrans_extra = (g_border_inner + pane_width)*col
        ytrans_initial = g_border_outer + pane_height/2.0
        ytrans_extra = (g_border_inner + pane_height)*row
        context.translate(
                xtrans_initial + xtrans_extra, ytrans_initial + ytrans_extra)
        # draw the tree into the context
        if bdrawbackground:
            bgcolor = g_color_background
        else:
            # Pretend that if there is not a background color
            # then the background color is white.
            bgcolor = (1.0, 1.0, 1.0)
        draw_single_tree(tree, context, v1, v2, bgcolor,
                bdrawvertices, bdrawlabels)
        # draw the pane label into the context
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
        # restore the position
        context.restore()
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

def _draw_bad_branches_revised(tree, context, v1, id_to_location):
    for node, child in tree.gen_directed_branches():
        if is_bad_edge_revised(node, child, v1):
            context.save()
            psrc = id_to_location[id(node)]
            pdst = id_to_location[id(child)]
            draw_wavy_line(
                    context,
                    psrc[0], psrc[1], pdst[0], pdst[1],
                    g_wavelength)
            context.restore()

def _draw_bad_branches(tree, context, v1, v2):
    for node, child in tree.gen_directed_branches():
        if is_bad_edge(node, child, v1, v2):
            context.save()
            psrc = tree._layout_to_display(node.location)
            pdst = tree._layout_to_display(child.location)
            #FIXME
            if False:
                context.set_source_rgb(*g_color_bad_edge)
                context.move_to(*psrc)
                context.line_to(*pdst)
                context.stroke()
            else:
                draw_wavy_line(
                        context,
                        psrc[0], psrc[1], pdst[0], pdst[1],
                        g_wavelength)
            context.restore()

def _draw_directed_branches_revised(tree, context, v1, id_to_location):
    for node, child in tree.gen_directed_branches():
        if is_bad_edge_revised(node, child, v1):
            continue
        context.save()
        context.set_source_rgb(*g_color_light)
        # Get the valuations and (x,y) points of each node.
        vsrc, vdst = v1[id(node)], v1[id(child)]
        psrc = id_to_location[id(node)]
        pdst = id_to_location[id(child)]
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

def _draw_vertex_ticks_revised(tree, context, v2, id_to_location):
    """
    @param v2: maps node id to valuation for zero crossing ticks
    """
    eps = 1e-8
    for node in tree.preorder():
        if abs(v2[id(node)]) > eps:
            continue
        npos = 0
        nneg = 0
        nzero = 0
        for n in node.get_neighbors():
            if v2[id(n)] < -eps:
                nneg += 1
            elif v2[id(n)] > eps:
                npos += 1
            else:
                nzero += 1
        if npos + nneg + nzero < 3:
            continue
        if not npos:
            continue
        if not nneg:
            continue
        x, y = id_to_location[id(node)]
        r = g_vertex_dot_radius_thick
        context.save()
        context.set_source_rgb(*g_color_dark)
        context.arc(x, y, r, 0, 2 * math.pi)
        context.fill()
        context.restore()

def _draw_edge_ticks_revised(tree, context, v2, id_to_location):
    for node, child in tree.gen_directed_branches():
        if is_bad_edge_revised(node, child, v2):
            continue
        context.save()
        # Get the valuations and (x,y) points of each node.
        vsrc, vdst = v2[id(node)], v2[id(child)]
        psrc = id_to_location[id(node)]
        pdst = id_to_location[id(child)]
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
            eps = 1e-8
            if abs(vsrc) > eps and abs(vdst) > eps:
                context.set_line_width(g_line_width_thin)
                context.move_to(barbx1, barby1)
                context.line_to(barbx2, barby2)
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
            eps = 1e-8
            if abs(vsrc) < eps or abs(vdst) < eps:
                #FIXME
                if False:
                    draw_wavy_line(
                            context,
                            barbx1, barby1, barbx2, barby2,
                            g_wavelength)
            else:
                context.set_line_width(g_line_width_thin)
                context.move_to(barbx1, barby1)
                context.line_to(barbx2, barby2)
                context.stroke()
        context.restore()

def _draw_labels_revised(tree, context, id_to_location):
    for node in tree.preorder():
        label = node.get_name()
        if label:
            # get the parameters for the label
            theta = get_free_angle_revised(tree, node, id_to_location)
            x, y = id_to_location[id(node)]
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

def _draw_labels(tree, context):
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

def get_free_angle_revised(tree, node, id_to_location):
    """
    @param tree: something like a SpatialTree
    @param node: something like a SpatialNode
    @return: get the angle from the node to some free space
    """
    # FIXME this is implemented quadratic but it could be implemented linear.
    # get the list of angles away from the node of interest
    angles = []
    origin = id_to_location[id(node)]
    for a, b in tree.gen_bidirected_branches():
        if a is node:
            target = id_to_location[id(b)]
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

def draw_single_tree_revised(
        tree, context, v1, v2, bgcolor, bdrawlabels, id_to_location):
    """
    This is most likely called only from inside the module.
    @param tree: a fitted SpatialTree
    @param context: a cairo context with origin at tree center
    @param v1: maps node id to valuation for line thickness
    @param v2: maps node id to valuation for zero crossing ticks
    @param bgcolor: background color
    @param bdrawlabels: flag to draw labels
    """
    _draw_bad_branches_revised(tree, context, v1, id_to_location)
    _draw_directed_branches_revised(tree, context, v1, id_to_location)
    if v2:
        _draw_edge_ticks_revised(tree, context, v2, id_to_location)
        _draw_vertex_ticks_revised(tree, context, v2, id_to_location)
    if bdrawlabels:
        _draw_labels_revised(tree, context, id_to_location)

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
                #r = g_vertex_dot_radius_thick
                #fgcolor = g_color_bad_edge
                r = g_vertex_dot_radius_thin
                fgcolor = g_color_dark
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
        _draw_labels(tree, context)

class TikzContext:
    def __init__(self):
        self.depth = 0
        self.lines = []
        self.finished = False
        self.add_line('\\begin{tikzpicture}[x=1cm,y=1cm,yscale=-1]')
        self.depth += 1
        # add the wave style
        style = '[snake=snake,color=gray,line after snake=0mm]'
        self.add_line('\\tikzstyle mywbasic='+style)
        style = '[style=mywbasic,segment amplitude=0.2mm,segment length=0.8mm]'
        self.add_line('\\tikzstyle mywave='+style)
        # add the thick and thin line styles
        self.add_line('\\tikzstyle mythick=[color=gray,line width=0.06cm]')
        self.add_line('\\tikzstyle mythin=[color=gray]')
    def add_line(self, line):
        if self.finished:
            raise ValueError('tried to add a line to a finished tikz context')
        self.lines.append((' ' * self.depth) + line)
    def begin_pane(self, xshift, yshift):
        self.add_line(
                '\\begin{scope}[xshift=%.4fcm,yshift=%.4fcm]' % (
                    xshift, yshift))
        self.depth += 1
    def end_pane(self):
        self.depth -= 1
        self.add_line('\\end{scope}')
    def finish(self):
        if not self.finished:
            self.depth -= 1
            self.add_line('\\end{tikzpicture}')
        self.finished = True
    def get_text(self):
        if not self.finished:
            raise ValueError('unfinished')
        return '\n'.join(self.lines)
    def draw_wavy_line(self, x1, y1, x2, y2):
        self.add_line(
                '\\draw[style=mywave] (%.4f,%.4f) -- (%.4f,%.4f);' % (
                    x1, y1, x2, y2))
    def draw_thin_light_line(self, x1, y1, x2, y2):
        self.add_line(
                '\\draw[style=mythin] (%.4f,%.4f) -- (%.4f,%.4f);' % (
                    x1, y1, x2, y2))
    def draw_thick_light_line(self, x1, y1, x2, y2):
        self.add_line(
                '\\draw[style=mythick] (%.4f,%.4f) -- (%.4f,%.4f);' % (
                    x1, y1, x2, y2))
    def draw_dark_line(self, x1, y1, x2, y2):
        self.add_line(
                '\\draw (%.4f,%.4f) -- (%.4f,%.4f);' % (
                    x1, y1, x2, y2))
    def draw_dark_dot(self, x, y):
        self.add_line(
                '\\fill (%.4f,%.4f) circle (0.04cm);' % (x, y))

def _draw_labels_tikz(tree, context, id_to_location):
    """
    Use degree anchors for label placement.
    """
    for node in tree.preorder():
        label = node.get_name()
        if label:
            # get the parameters for the label
            theta = get_free_angle_revised(tree, node, id_to_location)
            x, y = id_to_location[id(node)]
            # draw the text relative to the location
            theta = math.atan2(math.sin(theta), -math.cos(theta))
            float_degree = ((theta % (2 * math.pi)) * 360) / (2 * math.pi)
            #float_degree = (theta * 360) / (2 * math.pi)
            degree = int(math.floor(float_degree))
            style = 'font=\\tiny,anchor=%s,inner sep=1pt' % degree
            context.add_line(
                    '\\node[%s] at (%.4f,%.4f) {%s};' % (
                        style, x, y, label))

def _draw_vertex_ticks_tikz(tree, context, v2, id_to_location):
    """
    @param v2: maps node id to valuation for zero crossing ticks
    """
    eps = 1e-8
    for node in tree.preorder():
        if abs(v2[id(node)]) > eps:
            continue
        npos = 0
        nneg = 0
        nzero = 0
        for n in node.get_neighbors():
            if v2[id(n)] < -eps:
                nneg += 1
            elif v2[id(n)] > eps:
                npos += 1
            else:
                nzero += 1
        if npos + nneg + nzero < 3:
            continue
        if not npos:
            continue
        if not nneg:
            continue
        x, y = id_to_location[id(node)]
        context.draw_dark_dot(x, y)

def _draw_edge_ticks_tikz(tree, context, v2, id_to_location):
    barb_radius_cm = 0.06
    for node, child in tree.gen_directed_branches():
        if is_bad_edge_revised(node, child, v2):
            continue
        # Get the valuations and (x,y) points of each node.
        vsrc, vdst = v2[id(node)], v2[id(child)]
        psrc = id_to_location[id(node)]
        pdst = id_to_location[id(child)]
        if vsrc * vdst < 0:
            # find the crossing point
            t = -vsrc / (vdst - vsrc)
            pmid = [psrc[i] + t*(pdst[i] - psrc[i]) for i in range(2)]
            theta = math.atan2(pdst[1]-psrc[1], pdst[0]-psrc[0])
            barbx1 = pmid[0] + barb_radius_cm * math.cos(theta + math.pi/2)
            barby1 = pmid[1] + barb_radius_cm * math.sin(theta + math.pi/2)
            barbx2 = pmid[0] + barb_radius_cm * math.cos(theta - math.pi/2)
            barby2 = pmid[1] + barb_radius_cm * math.sin(theta - math.pi/2)
            # set line thickness for barb
            eps = 1e-8
            if abs(vsrc) > eps and abs(vdst) > eps:
                context.draw_dark_line(barbx1, barby1, barbx2, barby2)

def _draw_directed_branches_tikz(tree, context, v1, id_to_location):
    for node, child in tree.gen_directed_branches():
        if is_bad_edge_revised(node, child, v1):
            continue
        # Get the valuations and (x,y) points of each node.
        vsrc, vdst = v1[id(node)], v1[id(child)]
        psrc = id_to_location[id(node)]
        pdst = id_to_location[id(child)]
        if vsrc < 0 and vdst < 0:
            context.draw_thin_light_line(
                    psrc[0], psrc[1], pdst[0], pdst[1])
        elif vsrc > 0 and vdst > 0:
            context.draw_thick_light_line(
                    psrc[0], psrc[1], pdst[0], pdst[1])
        else:
            # find the crossing point
            t = -vsrc / (vdst - vsrc)
            pmid = [psrc[i] + t*(pdst[i] - psrc[i]) for i in range(2)]
            # set line thickness for source
            if vsrc < 0:
                context.draw_thin_light_line(
                        psrc[0], psrc[1], pmid[0], pmid[1])
            else:
                context.draw_thick_light_line(
                        psrc[0], psrc[1], pmid[0], pmid[1])
            # set line thickness for destination
            if vdst < 0:
                context.draw_thin_light_line(
                        pmid[0], pmid[1], pdst[0], pdst[1])
            else:
                context.draw_thick_light_line(
                        pmid[0], pmid[1], pdst[0], pdst[1])

def _draw_bad_branches_tikz(tree, context, v1, id_to_location):
    for node, child in tree.gen_directed_branches():
        if is_bad_edge_revised(node, child, v1):
            psrc = id_to_location[id(node)]
            pdst = id_to_location[id(child)]
            context.draw_wavy_line(psrc[0], psrc[1], pdst[0], pdst[1])

def draw_single_tree_tikz(tree, context, v1, v2,
        flag_draw_labels, id_to_location):
    """
    This is most likely called only from inside the module.
    @param tree: a fitted SpatialTree
    @param context: a tikz context with origin at tree center
    @param v1: maps node id to valuation for line thickness
    @param v2: maps node id to valuation for zero crossing ticks
    """
    _draw_bad_branches_tikz(tree, context, v1, id_to_location)
    _draw_directed_branches_tikz(tree, context, v1, id_to_location)
    if v2:
        _draw_edge_ticks_tikz(tree, context, v2, id_to_location)
        _draw_vertex_ticks_tikz(tree, context, v2, id_to_location)
    if flag_draw_labels:
        _draw_labels_tikz(tree, context, id_to_location)

def get_forest_image_tikz_old(
            tree, max_size, vs, inner_margin,
            reflect_trees, flag_draw_labels):
    """
    Get the image of the tree.
    Use a manual grid layout.
    This could be called from outside the module.
    The rectangle size is determined automatically
    so as to maximize the scaling factors of the trees in the panels.
    @param tree: something like a SpatialTree
    @param max_size: (max_width, max_height) in centimeters
    @param vs: sequence of maps from node id to valuation
    @param inner_margin: inner margin in centimeters
    @param reflect_trees: a flag for tweaking the tree layout
    @return: tikz text
    """
    outer_margin = 0
    npairs = len(vs)
    if npairs < 1:
        raise ValueError('not enough valuation maps')
    # get all minimum rectangle sizes
    rect_sizes = layout.get_rect_sizes(npairs)
    # get (scaling_factor, w, h) triples
    triples = []
    max_width, max_height = max_size
    for ncols, nrows in rect_sizes:
        max_pane_width = (
                max_width - 2*outer_margin - inner_margin*(ncols-1))/ncols
        max_pane_height = (
                max_height - 2*outer_margin - inner_margin*(nrows-1))/nrows
        # require a minimum size
        min_pane_centimeters = 1e-4
        if max_pane_width < min_pane_centimeters:
            continue
        elif max_pane_height < min_pane_centimeters:
            continue
        # append a triple after finding the scaling factor
        max_pane_size = (max_pane_width, max_pane_height)
        tree.fit(max_pane_size)
        triples.append((tree.scale, ncols, nrows))
    # if no triples were found then we fail
    if not triples:
        raise ValueError('not enough room')
    # get the nrows and ncols corresponding to the best scaling factor
    max_scale, ncols, nrows = max(triples)
    # get the position of each pane
    row_col_pairs = layout.min_rect_to_row_major(ncols, nrows, npairs)
    # re-fit the tree using the best size
    max_pane_width = (
            max_width - 2*outer_margin - inner_margin*(ncols-1))/ncols
    max_pane_height = (
            max_height - 2*outer_margin - inner_margin*(nrows-1))/nrows
    max_pane_size = (max_pane_width, max_pane_height)
    tree.fit(max_pane_size)
    # get the map from id to location for the final tree layout
    nodes = list(tree.preorder())
    id_to_location = dict(
            (id(n), tree._layout_to_display(n.location)) for n in nodes)
    if reflect_trees:
        id_to_location = dict(
                (v, (-x, y)) for v, (x,y) in id_to_location.items())
    # get the width and height of the tree image
    xmin, ymin, xmax, ymax = locations_to_extents(id_to_location.values())
    pane_width = xmax - xmin
    pane_height = ymax - ymin
    width = 2*outer_margin + (ncols-1)*inner_margin + ncols*pane_width
    height = 2*outer_margin + (nrows-1)*inner_margin + nrows*pane_height
    # draw the trees
    context = TikzContext()
    for i, ((row, col), (v1, v2)) in enumerate(
            zip(row_col_pairs, iterutils.pairwise(vs+[None]))):
        # center the context on the correct pane
        xtrans_initial = outer_margin + pane_width/2.0
        xtrans_extra = (inner_margin + pane_width)*col
        ytrans_initial = outer_margin + pane_height/2.0
        ytrans_extra = (inner_margin + pane_height)*row
        xshift = xtrans_initial + xtrans_extra
        yshift = ytrans_initial + ytrans_extra
        context.begin_pane(xshift, yshift)
        # draw the tree into the context
        draw_single_tree_tikz(
                tree, context, v1, v2, flag_draw_labels, id_to_location)
        # draw the pane label into the context
        if i < len(string.uppercase):
            pane_label = str(i+1)
        else:
            pane_label = '?'
        xtarget = -pane_width/2
        ytarget = -pane_height/2
        style = 'anchor=north west,inner sep=1pt'
        context.add_line(
            '\\node[%s] at (%.4f,%.4f) {%s};' % (
                style, xtarget, ytarget, pane_label))
        context.end_pane()
    context.finish()
    return context.get_text()

def get_forest_image_tikz(
            tree, max_size, vs, inner_margin,
            reflect_trees, flag_draw_labels):
    """
    Get the image of the tree.
    Attempt to use the builtin TikZ matrix layout.
    This could be called from outside the module.
    The rectangle size is determined automatically
    so as to maximize the scaling factors of the trees in the panels.
    @param tree: something like a SpatialTree
    @param max_size: (max_width, max_height) in centimeters
    @param vs: sequence of maps from node id to valuation
    @param inner_margin: inner margin in centimeters
    @param reflect_trees: a flag for tweaking the tree layout
    @return: tikz text
    """
    outer_margin = 0
    npairs = len(vs)
    if npairs < 1:
        raise ValueError('not enough valuation maps')
    # get all minimum rectangle sizes
    rect_sizes = layout.get_rect_sizes(npairs)
    # get (scaling_factor, w, h) triples
    triples = []
    max_width, max_height = max_size
    for ncols, nrows in rect_sizes:
        max_pane_width = (
                max_width - 2*outer_margin - inner_margin*(ncols-1))/ncols
        max_pane_height = (
                max_height - 2*outer_margin - inner_margin*(nrows-1))/nrows
        # require a minimum size
        min_pane_centimeters = 1e-4
        if max_pane_width < min_pane_centimeters:
            continue
        elif max_pane_height < min_pane_centimeters:
            continue
        # append a triple after finding the scaling factor
        max_pane_size = (max_pane_width, max_pane_height)
        tree.fit(max_pane_size)
        triples.append((tree.scale, ncols, nrows))
    # if no triples were found then we fail
    if not triples:
        raise ValueError('not enough room')
    # get the nrows and ncols corresponding to the best scaling factor
    max_scale, ncols, nrows = max(triples)
    # get the position of each pane
    row_col_pairs = layout.min_rect_to_row_major(ncols, nrows, npairs)
    # re-fit the tree using the best size
    max_pane_width = (
            max_width - 2*outer_margin - inner_margin*(ncols-1))/ncols
    max_pane_height = (
            max_height - 2*outer_margin - inner_margin*(nrows-1))/nrows
    max_pane_size = (max_pane_width, max_pane_height)
    tree.fit(max_pane_size)
    # get the map from id to location for the final tree layout
    nodes = list(tree.preorder())
    id_to_location = dict(
            (id(n), tree._layout_to_display(n.location)) for n in nodes)
    if reflect_trees:
        id_to_location = dict(
                (v, (-x, y)) for v, (x,y) in id_to_location.items())
    # get the width and height of the tree image
    xmin, ymin, xmax, ymax = locations_to_extents(id_to_location.values())
    pane_width = xmax - xmin
    pane_height = ymax - ymin
    width = 2*outer_margin + (ncols-1)*inner_margin + ncols*pane_width
    height = 2*outer_margin + (nrows-1)*inner_margin + nrows*pane_height
    # draw the trees
    context = TikzContext()
    style = 'column sep=%.4f,row sep=%.4f' % (
            inner_margin, inner_margin)
    context.add_line('\\matrix[%s] {' % style)
    context.depth += 1
    for i, ((row, col), (v1, v2)) in enumerate(
            zip(row_col_pairs, iterutils.pairwise(vs+[None]))):
        # draw the tree into the context
        context.add_line('\\begin{scope}[yscale=-1]')
        context.depth += 1
        draw_single_tree_tikz(
                tree, context, v1, v2, flag_draw_labels, id_to_location)
        # draw the pane label into the context
        if i < len(string.uppercase):
            pane_label = str(i+1)
        else:
            pane_label = '?'
        xtarget = -pane_width/2
        ytarget = -pane_height/2
        style = 'anchor=north west,inner sep=1pt'
        context.add_line(
            '\\node[%s] at (%.4f,%.4f) {%s};' % (
                style, xtarget, ytarget, pane_label))
        context.depth -= 1
        context.add_line('\\end{scope}')
        # add a row break or a column break
        nblanks = nrows * ncols - len(row_col_pairs)
        if i == len(row_col_pairs) - 1:
            line = ''
            if nblanks:
                line += ' '.join(['&']*nblanks)
            line += '\\\\'
            context.add_line(line)
        elif col == ncols-1:
            context.add_line('\\\\')
        else:
            context.add_line('&')
    context.depth -= 1
    context.add_line('};')
    context.finish()
    return context.get_text()

