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


from StringIO import StringIO
import math
import random
import itertools
import time

import numpy as np
from matplotlib.delaunay.triangulate import Triangulation

from SnippetUtil import HandlingError
import SnippetUtil
import Form
import FormOut
import CompGeom

#FIXME clarify output format

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
                10.0, low_exclusive=0.0),
            Form.Float('ann_sigma', 'stddev of annulus',
                1.2, low_exclusive=0.0),
            Form.Float('disc_sigma', 'stddev of disc',
                1.0, low_exclusive=0.0),
            Form.Float('conn_radius', 'remove edges longer than this',
                5.0, low_exclusive=0.0),
            Form.Integer('conn_sparse', 'keep this many disc-annulus edges',
                2, low=0),
            Form.Integer('disc_npoints', 'number of points in the disc',
                50, low=1, high=500),
            Form.Integer('ann_npoints', 'number of points in the annulus',
                150, low=1, high=500)]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # start writing the response
    out = StringIO()
    # define the points; each is a numpy array of length 2
    G = []
    P = []
    for i in range(fs.disc_npoints):
        P.append(sample_point_on_disc(fs.disc_sigma))
        G.append(0)
    for i in range(fs.ann_npoints):
        P.append(sample_point_on_annulus(fs.ann_radius, fs.ann_sigma))
        G.append(1)
    # use delaunay triangulation to force planarity
    x_list, y_list = zip(*[p.tolist() for p in P])
    tri = Triangulation(x_list, y_list)
    edges = set(tuple(sorted(edge)) for edge in tri.edge_db)
    # define edge distances
    e_to_d = dict(((i, j), np.linalg.norm(P[j]-P[i])) for i, j in edges)
    # get the edges joining the disc and annulus groups
    joining_edges = set((i,j) for i, j in edges if G[i] != G[j])
    dae_pairs = [(e_to_d[e], e) for e in joining_edges]
    sorted_dae_pairs = list(sorted(dae_pairs))
    sorted_joining_edges = zip(*sorted_dae_pairs)[1]
    short_joining_edges = set(sorted_joining_edges[:fs.conn_sparse])
    # get edges which are longer than the connection radius
    overlong_edges = set(e for e in edges if e_to_d[e] > fs.conn_radius)
    # define the final set of edges
    edges = (edges - overlong_edges - joining_edges) | short_joining_edges
    # write some extra info
    # write the points
    print >> out, 'POINTS'
    for i, p in enumerate(P):
        print >> out, '\t'.join(str(x) for x in [i, p[0], p[1]])
    # write the edges
    print >> out, 'EDGES'
    for i, j in edges:
        print >> out, '\t'.join(str(x) for x in [i, j])
    # write the response
    return [('Content-Type', 'text/plain')], out.getvalue().strip()
