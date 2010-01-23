"""Sample a disk-annulus mixture of points and define some edges.

The output of this snippet is supposed to be used for spectral partitioning.
The output format is as follows.
POINTS
index_0 x0 y0
index_1 x1 y1
index_n xn yn
EDGES
index_a1 index_b1
index_a2 index_b2
index_ak index_bk
"""


import StringIO
import math
import random
import itertools
import time

import numpy as np
from matplotlib.delaunay.triangulate import Triangulation

from SnippetUtil import HandlingError
import SnippetUtil
import Form
import CompGeom


def sample_point_on_circle(radius):
    """
    The circle is centered at the origin and has the given radius.
    @param radius: the radius of the circle
    @return: a uniformly random point on the circle
    """
    theta = random.random() * 2 * math.pi
    x = radius * math.cos(theta)
    y = radius * math.sin(theta)
    return np.array([x, y])

def sample_point_on_annulus(radius, sigma):
    """
    The annulus is centered at the origin.
    This is not really an annulus but rather a fuzzy circle.
    @param radius: the radius of center of the annulus
    @param sigma: the fuzziness of the circle
    @return: a sampled point on the annulus
    """
    p = sample_point_on_circle(radius)
    dx = random.gauss(0, sigma)
    dy = random.gauss(0, sigma)
    return p + np.array([dx, dy])

def sample_point_on_disc(sigma):
    """
    This is not really a disc but rather a fuzzy point.
    @param sigma: the fuzziness of the point
    @return: a sampled point
    """
    x = random.gauss(0, sigma)
    y = random.gauss(0, sigma)
    return np.array([x, y])

def get_intersecting_edges(points, edges):
    """
    Edges that share an endpoint do not count as conflicting.
    @param points: a list of numpy arrays each of length 2
    @param edges: a list of point index pairs
    @return: the set of edges that intersect at least one other edge
    """
    conflicts = set()
    sedgewick_points = [CompGeom.Point(p.tolist()) for p in points]
    for ea, eb in itertools.combinations(edges, 2):
        # only check intersection when each endpoint is unique
        if len(set([ea[0], ea[1], eb[0], eb[1]])) == 4:
            a = sedgewick_points[ea[0]]
            b = sedgewick_points[ea[1]]
            c = sedgewick_points[eb[0]]
            d = sedgewick_points[eb[1]]
            if CompGeom.line_segments_intersect(a, b, c, d):
                conflicts.add(ea)
                conflicts.add(eb)
    return conflicts

def get_form():
    """
    @return: a list of form objects
    """
    # define the list of form objects
    form_objects = [
            Form.Float('ann_radius', 'radius of annulus center',
                3.0, low_exclusive=0.0),
            Form.Float('ann_sigma', 'stddev of annulus',
                1.0, low_exclusive=0.0),
            Form.Float('disc_sigma', 'stddev of disc',
                1.0, low_exclusive=0.0),
            Form.Float('conn_radius', 'connection radius',
                0.2, low_exclusive=0.0),
            Form.Integer('disc_npoints', 'number of points in the disc',
                20, low=1, high=500),
            Form.Integer('ann_npoints', 'number of points in the annulus',
                20, low=1, high=500),
            Form.CheckGroup('options', 'graph options', [
                Form.CheckItem('force_planar', 'force planar graph', True)])]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # start writing the response
    out = StringIO.StringIO()
    # define the points; each is a numpy array of length 2
    points = []
    for i in range(fs.disc_npoints):
        points.append(sample_point_on_disc(fs.disc_sigma))
    for i in range(fs.ann_npoints):
        points.append(sample_point_on_annulus(fs.ann_radius, fs.ann_sigma))
    # create some short edges; each is a pair of point indices
    short_edges = set()
    for i, pa in enumerate(points):
        for j, pb in enumerate(points):
            if i < j:
                if np.linalg.norm(pb-pa) < fs.conn_radius:
                    edge = (i, j)
                    short_edges.add(edge)
    # use delaunay triangulation to force planarity
    if fs.force_planar:
        x_list, y_list = zip(*[p.tolist() for p in points])
        tri = Triangulation(x_list, y_list)
        delaunay_edges = set(tuple(sorted(edge)) for edge in tri.edge_db)
        edges = short_edges & delaunay_edges
    else:
        edges = short_edges
    # write some extra info
    # write the points
    print >> out, 'POINTS'
    for i, p in enumerate(points):
        print >> out, '\t'.join(str(x) for x in [i, p[0], p[1]])
    # write the edges
    print >> out, 'EDGES'
    for i, j in edges:
        print >> out, '\t'.join(str(x) for x in [i, j])
    # write the response
    return [('Content-Type', 'text/plain')], out.getvalue().strip()
