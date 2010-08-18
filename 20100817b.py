"""Sample a minimum spanning tree on a plane.
"""

import random
import math

import numpy as np
import cairo
from matplotlib.delaunay.triangulate import Triangulation

from SnippetUtil import HandlingError
import Form
import FormOut
import CompGeom
import CairoUtil
import const

g_tree_string = const.read('20100730g').rstrip()

g_africa_poly = [
        (686., 274.),
        (653., 368.),
        (576., 462.),
        (589., 572.),
        (533., 619.),
        (533., 668.),
        (447., 768.),
        (380., 782.),
        (306., 594.),
        (326., 527.),
        (281., 418.),
        (282., 363.),
        (252., 363.),
        (234., 339.),
        (124., 363.),
        (36., 255.),
        (40., 173.),
        (157., 18.),
        (296., 3.),
        (298., 41.),
        (377., 78.),
        (398., 48.),
        (512., 68.),
        (617., 289.)]


def get_form():
    """
    @return: the body of a form
    """
    # define the form objects
    form_objects = [
            Form.Integer('npoints', 'sample this many points', 200),
            Form.Integer('total_width', 'total image width',
                640, low=3, high=2000),
            Form.Integer('total_height', 'total image height',
                480, low=3, high=2000),
            Form.Integer('border', 'image border size',
                10, low=0, high=2000),
            Form.ImageFormat(),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.Image('demo')

def sample_with_rejection(npoints, poly, limit):
    """
    @param npoints: sample this many points with rejection
    @param poly: reject points outside the polygon
    @param limit: a sanity limit for rejection sampling
    @return: the sampled points
    """
    xpoly, ypoly = zip(*poly)
    xlow, xhigh = min(xpoly), max(xpoly)
    ylow, yhigh = min(ypoly), max(ypoly)
    points = []
    for i in range(limit):
        x = random.uniform(xlow, xhigh)
        y = random.uniform(ylow, yhigh)
        if CompGeom.point_in_poly(x, y, poly):
            points.append((x,y))
        if len(points) == npoints:
            break
    return points

def gen_noncrossing_edges(query_edges, boundary_edges, allpoints):
    """
    Yield query edges which do not intersect boundary edges.
    @param query_edges: point index pairs
    @param boundary_edges: point index pairs
    @param allpoints: ordered points as pairs of coordinates
    """
    for query in query_edges:
        q0 = CompGeom.Point(allpoints[query[0]])
        q1 = CompGeom.Point(allpoints[query[1]])
        for boundary in boundary_edges:
            b0 = CompGeom.Point(allpoints[boundary[0]])
            b1 = CompGeom.Point(allpoints[boundary[1]])
            if CompGeom.line_segments_intersect(q0, q1, b0, b1):
                break
        else:
            yield query

def get_response_content(fs):
    limit = fs.npoints*100
    # get points defining the boundary of africa
    nafrica = len(g_africa_poly)
    africa_edges = [(i, (i+1)%nafrica) for i in range(nafrica)]
    # get some points and edges inside africa
    points = sample_with_rejection(fs.npoints, g_africa_poly, limit)
    x_list, y_list = zip(*points)
    tri = Triangulation(x_list, y_list)
    tri_edges = [(i+nafrica, j+nafrica) for i, j in tri.edge_db.tolist()]
    # get the whole list of points
    allpoints = g_africa_poly + points
    # refine the list of edges
    tri_edges = list(gen_noncrossing_edges(tri_edges, africa_edges, allpoints))
    alledges = africa_edges + tri_edges
    # get the width and height of the drawable area of the image
    width = fs.total_width - 2*fs.border
    height = fs.total_height - 2*fs.border
    if width < 1 or height < 1:
        msg = 'the image dimensions do not allow for enough drawable area'
        raise HandlingError(msg)
    # draw the image
    try:
        ext = Form.g_imageformat_to_ext[fs.imageformat]
        return get_image_string(allpoints, alledges,
                fs.total_width, fs.total_height, fs.border, ext)
    except CairoUtil.CairoUtilError, e:
        raise HandlingError(e)

def get_image_string(points, edges, t_width, t_height, border, image_format):
    """
    @param points: an ordered list of (x, y) pairs
    @param edges: a set of point index pairs
    @param t_width: the image width in pixels
    @param t_height: the image height in pixels
    @param border: the width and height of the image border in pixels
    @param image_format: the requested format of the image
    """
    width = t_width - 2*border
    height = t_height - 2*border
    assert width >= 1
    assert height >= 1
    # get the x and y coordinates of the points
    x_coords, y_coords = zip(*points)
    unscaled_width = max(x_coords) - min(x_coords)
    unscaled_height = max(y_coords) - min(y_coords)
    # rescale the points to fit in the drawable part of the image
    c_width = width / unscaled_width
    c_height = height / unscaled_height
    c = min(c_width, c_height)
    x_rescaled = [x*c for x in x_coords]
    y_rescaled = [y*c for y in y_coords]
    # translate the points so that their minimum coordinate is zero
    x_min = min(x_rescaled)
    y_min = min(y_rescaled)
    x_trans = [x-x_min for x in x_rescaled]
    y_trans = [y-y_min for y in y_rescaled]
    # translate the points to account for the border
    x_final = [x+border for x in x_trans]
    y_final = [y+border for y in y_trans]
    # create the surface
    cairo_helper = CairoUtil.CairoHelper(image_format)
    surface = cairo_helper.create_surface(t_width, t_height)
    context = cairo.Context(surface)
    # draw the background
    context.save()
    context.set_source_rgb(.9, .9, .9)
    context.paint()
    context.restore()
    # draw the points
    radius = 2.0
    for x, y in zip(x_final, y_final):
        # draw a filled circle
        context.save()
        context.arc(x, y, radius, 0, 2 * math.pi)
        context.close_path()
        context.fill()
        context.restore()
    # draw the edges
    for i, j in edges:
        context.save()
        context.move_to(x_final[i], y_final[i])
        context.line_to(x_final[j], y_final[j])
        context.close_path()
        context.stroke()
        context.restore()
    # get the image string
    return cairo_helper.get_image_string()
