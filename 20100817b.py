"""Sample a minimum spanning tree on a plane.
"""

import random
import math
import itertools

import numpy as np
import cairo
from matplotlib.delaunay.triangulate import Triangulation
import matplotlib.pyplot as plt
import matplotlib
import argparse

from SnippetUtil import HandlingError
import Form
import FormOut
import CompGeom
import CairoUtil
import EigUtil
import MST
import const
import Euclid

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
            Form.Integer('seed', 'rng seed or 0 for no fixed seed', 0),
            Form.Integer('npoints', 'sample this many points', 200),
            Form.Integer('axis', 'color according to this MDS axis', 0),
            Form.Integer('total_width', 'total image width',
                640, low=3, high=2000),
            Form.Integer('total_height', 'total image height',
                480, low=3, high=2000),
            Form.Integer('border', 'image border size',
                10, low=0, high=2000),
            Form.ImageFormat()]
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

def get_mst(query_edges, allpoints):
    """
    Return a sublist of the query edges which forms a minimum spanning tree.
    @param query_edges: point index pairs
    @param allpoints: ordered points as pairs of coordinates
    @return: minimum spanning tree edges
    """
    # get the point indices occurring in at least one query edge
    V = list(sorted(set(itertools.chain.from_iterable(query_edges))))
    # get the weighted edges
    E = []
    for edge in query_edges:
        ax, ay = allpoints[edge[0]]
        bx, by = allpoints[edge[1]]
        weight = math.hypot(bx-ax, by-ay)
        E.append((weight, edge[0], edge[1]))
    # get the minimum spanning tree weighted edges
    E_mst = MST.kruskal(V, E)
    # get the minimum spanning tree edges
    return [(va, vb) for w, va, vb in E_mst]

def get_color(x):
    """
    @param x: between -1 and 1
    @return: an rgb triple in [0,1]
    """
    r = 1.0 if x > 0 else 1.0 - abs(x)
    g = 1.0 if x < 0 else 1.0 - abs(x)
    b = 1.0 - abs(x)
    return r, g, b

def get_response_content(fs):
    # use a fixed seed if requested
    if fs.seed:
        random.seed(fs.seed)
    # define the max number of rejection iterations
    limit = fs.npoints*100
    # validate input
    if fs.axis < 0:
        raise ValueError('the mds axis must be nonnegative')
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
    tri_edges = get_mst(tri_edges, allpoints)
    alledges = africa_edges + tri_edges
    # make the graph laplacian
    A = np.zeros((len(points), len(points)))
    for ia, ib in tri_edges:
        xa, ya = allpoints[ia]
        xb, yb = allpoints[ib]
        d = math.hypot(xb-xa, yb-ya)
        A[ia-nafrica, ib-nafrica] = 1/d
        A[ib-nafrica, ia-nafrica] = 1/d
    L = Euclid.adjacency_to_laplacian(A)
    ws, vs = EigUtil.eigh(np.linalg.pinv(L))
    if fs.axis >= len(ws):
        raise ValueError('choose a smaller mds axis')
    v = vs[fs.axis]
    # get the color and sizes for the points
    v /= max(np.abs(v))
    colors = [(0,0,0)]*nafrica + [get_color(x) for x in v]
    radii = [2]*nafrica + [5 for p in points]
    # get the width and height of the drawable area of the image
    width = fs.total_width - 2*fs.border
    height = fs.total_height - 2*fs.border
    if width < 1 or height < 1:
        msg = 'the image dimensions do not allow for enough drawable area'
        raise HandlingError(msg)
    # draw the image
    ext = Form.g_imageformat_to_ext[fs.imageformat]
    try:
        helper = ImgHelper(allpoints, alledges,
                fs.total_width, fs.total_height, fs.border)
        return helper.get_image_string(colors, radii, ext)
    except CairoUtil.CairoUtilError as e:
        raise HandlingError(e)


class ImgHelper:
    """
    This used to be a single function but there is too much going on in it.
    """

    def __init__(self, points, edges, t_width, t_height, border):
        """
        @param points: an ordered list of (x, y) pairs
        @param edges: a set of point index pairs
        @param t_width: the image width in pixels
        @param t_height: the image height in pixels
        @param border: the width and height of the image border in pixels
        """
        self.edges = edges
        self.t_width = t_width
        self.t_height = t_height
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
        self.x_final = [x+border for x in x_trans]
        self.y_final = [y+border for y in y_trans]

    def draw_contour_plot(self, valuations, nafrica):
        """
        Note that this is not useful for the web.
        Also it is not useful for calling over ssh
        because it pops up a picture.
        @param valuations: valuations conformant to the spanning tree points
        @param nafrica: the number of points on the continent boundary
        """
        # y axis is flipped in matplotlib relative to cairo
        x = self.x_final
        y = [self.t_height - y for y in self.y_final]
        # get the relevant points
        mst_x = x[nafrica:]
        mst_y = y[nafrica:]
        # do a grid interpolation
        xvalues = np.arange(min(mst_x), max(mst_x), 1.0)
        yvalues = np.arange(min(mst_y), max(mst_y), 1.0)
        xi, yi = np.meshgrid(xvalues, yvalues)
        tri = Triangulation(mst_x, mst_y)
        interp = tri.nn_interpolator(valuations)
        zi = interp(xi, yi)
        # define the africa clip path
        # http://matplotlib.sourceforge.net/examples/api/compound_path.html
        continent_poly = np.array(zip(x[:nafrica], y[:nafrica]) + [(0,0)])
        continent_codes = ([matplotlib.path.Path.MOVETO] +
                [matplotlib.path.Path.LINETO]*(nafrica-1) +
                [matplotlib.path.Path.CLOSEPOLY])
        continent_path = matplotlib.path.Path(continent_poly, continent_codes)
        # draw the plot
        plt.figure()
        cs = plt.contour(xi, yi, zi, clip_path=continent_path)
        plt.fill(x[:nafrica], y[:nafrica],
                fill=False)
        plt.axes().set_aspect('equal')
        plt.show()

    def get_image_string(self, colors, radii, image_format):
        """
        @param colors: point colors as rgb triples in [0,1]
        @param radii: point radii should be about 2.0
        @param image_format: the requested format of the image
        @return: the image string
        """
        # create the surface
        cairo_helper = CairoUtil.CairoHelper(image_format)
        surface = cairo_helper.create_surface(self.t_width, self.t_height)
        context = cairo.Context(surface)
        # draw the background
        context.save()
        context.set_source_rgb(.9, .9, .9)
        context.paint()
        context.restore()
        # draw the edges
        for i, j in self.edges:
            context.save()
            context.move_to(self.x_final[i], self.y_final[i])
            context.line_to(self.x_final[j], self.y_final[j])
            context.close_path()
            context.stroke()
            context.restore()
        # draw the points
        for r, c, x, y in zip(radii, colors, self.x_final, self.y_final):
            # draw a filled circle
            context.save()
            context.set_source_rgb(*c)
            context.arc(x, y, r, 0, 2 * math.pi)
            context.close_path()
            context.fill()
            context.restore()
        # get the image string
        return cairo_helper.get_image_string()

def main(fs):
    # use a fixed seed if requested
    if fs.seed:
        random.seed(fs.seed)
    # define the max number of rejection iterations
    limit = fs.npoints*100
    # validate input
    if fs.axis < 0:
        raise ValueError('the mds axis must be nonnegative')
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
    tri_edges = get_mst(tri_edges, allpoints)
    alledges = africa_edges + tri_edges
    # make the graph laplacian
    A = np.zeros((len(points), len(points)))
    for ia, ib in tri_edges:
        xa, ya = allpoints[ia]
        xb, yb = allpoints[ib]
        d = math.hypot(xb-xa, yb-ya)
        A[ia-nafrica, ib-nafrica] = 1/d
        A[ib-nafrica, ia-nafrica] = 1/d
    L = Euclid.adjacency_to_laplacian(A)
    ws, vs = EigUtil.eigh(np.linalg.pinv(L))
    if fs.axis >= len(ws):
        raise ValueError('choose a smaller mds axis')
    v = vs[fs.axis]
    # get the color and sizes for the points
    v /= max(np.abs(v))
    # draw the picture
    helper = ImgHelper(allpoints, alledges,
            fs.total_width, fs.total_height, fs.border)
    helper.draw_contour_plot(v, nafrica)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--seed', default=0, type=int,
            help='random number seed or 0')
    parser.add_argument('--axis', default=0, type=int,
            help='color using this axis')
    parser.add_argument('--npoints', default=100, type=int,
            help='sample this many points')
    parser.add_argument('--total_width', default=720, type=int,
            help='image width')
    parser.add_argument('--total_height', default=720, type=int,
            help='image height')
    parser.add_argument('--border', default=10, type=int,
            help='image border')
    main(parser.parse_args())
