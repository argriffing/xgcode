#!/usr/bin/env python

"""Implement some computational geometry algorithms.

Many of these functions are inspired by Algorithms in C++ by Robert Sedgewick.
"""

import math
import random
import unittest
import copy

class Point:
    def __init__(self, pair=None):
        if pair:
            self.x, self.y = pair
    def to_tuple(self):
        return (self.x, self.y)
    def __repr__(self):
        return 'Point(%s)' % str(self.to_tuple())


class ConvexPolygonError(Exception):
    pass


class ConvexPolygon:
    """
    This is a glorified list of points.
    Following the vertices should give a counter clockwise path.
    """
    def __init__(self, points=()):
        self.points = points
    def assert_valid(self):
        if len(self.points) < 3:
            raise ConvexPolygonError('too few points: %d' % len(self.points))
        extended = self.points + self.points[:2]
        triples = zip(extended[:-2], extended[1:-1], extended[2:])
        for a, b, c in triples:
            result = ccw(a, b, c)
            if result == 0:
                raise ConvexPolygonError('found collinear points')
            elif result == -1:
                raise ConvexPolygonError('found clockwise points')
    def gen_edges(self):
        """
        Yield ordered point pairs representing directed edges.
        """
        extended = self.points + self.points[:1]
        for pair in zip(extended[:-1], extended[1:]):
            yield pair
    def covers_point(self, point):
        """
        @param point: an object with x and y members
        @return: True when the point is in the convex polygon
        """
        # check vertex collision
        for vertex in self.points:
            if (point.x, point.y) == vertex.to_tuple():
                return False
        n = len(self.points)
        for b, c in self.gen_edges():
            result = ccw(point, b, c)
            if result == 0:
                raise ConvexPolygonError('found collinear points')
            elif result == -1:
                return False
        return True


class CaliperState:
    
    def __init__(self, convex_hull=None):
        if convex_hull:
            xmin_index = min((p.x, i) for i, p in enumerate(convex_hull))[1]
            ymin_index = min((p.y, i) for i, p in enumerate(convex_hull))[1]
            xmax_index = max((p.x, i) for i, p in enumerate(convex_hull))[1]
            ymax_index = max((p.y, i) for i, p in enumerate(convex_hull))[1]
            self.corner_indices = [ymin_index, xmax_index, ymax_index, xmin_index]
            self.rotation = 0
            self.convex_hull = convex_hull

    def clone(self):
        """
        Do a somewhat complex copy that is neither shallow nor deep.
        """
        new = CaliperState()
        # copy the corner indices and the rotation
        new.corner_indices = self.corner_indices[:]
        new.rotation = self.rotation
        # link to the convex hull
        new.convex_hull = self.convex_hull
        return new

    def gen_rectangle_vertices(self):
        """
        Yield ordered vertices of the bounding rectangle defined by the calipers.
        """
        for i, index in enumerate(self.corner_indices):
            p0 = self.convex_hull[index]
            p1 = self.convex_hull[self.corner_indices[(i+1)%4]]
            current_angle = self.rotation + i * (math.pi / 2)
            dy = p1.y - p0.y
            dx = p1.x - p0.x
            distance_to_next = math.hypot(dy, dx)
            angle_to_next = math.atan2(dy, dx)
            angle_difference = (angle_to_next - current_angle) % (2*math.pi)
            extension_distance = distance_to_next * math.cos(angle_difference)
            x = p0.x + extension_distance * math.cos(current_angle)
            y = p0.y + extension_distance * math.sin(current_angle)
            yield Point((x, y))

    def get_next_state(self):
        """
        @return: the caliper state after rotating to the next edge
        """
        npoints = len(self.convex_hull)
        rotation_corner_pairs = []
        for i, index in enumerate(self.corner_indices):
            p0 = self.convex_hull[index]
            p1 = self.convex_hull[(index+1)%npoints]
            current_angle = self.rotation + i * (math.pi / 2)
            next_angle = math.atan2(p1.y - p0.y, p1.x - p0.x)
            angle_difference = (next_angle - current_angle) % (2*math.pi)
            rotation_corner_pairs.append((angle_difference, i))
        best_rotation, best_corner = min(rotation_corner_pairs)
        next_state = self.clone()
        next_state.rotation = self.rotation + best_rotation
        next_state.corner_indices[best_corner] = (next_state.corner_indices[best_corner] + 1) % npoints
        return next_state

def line_segments_intersect(a, b, c, d):
    """
    Decide whether a pair of line segments intersects.
    @param a: the first endpoint of the first line segment
    @param b: the second endpoint of the first line segment
    @param c: the first endpoint of the second line segment
    @param d: the second endpoint of the second line segment
    @return: True if line segments intersect
    """
    if ccw(a,c,d) == ccw(b,c,d):
        return False
    if ccw(a,b,c) == ccw(a,b,d):
        return False
    return True

def ccw(p0, p1, p2):
    """
    @return: 0 when collinear, 1 when turning ccw, -1 when turning cw
    """
    dx1 = p1.x - p0.x
    dy1 = p1.y - p0.y
    dx2 = p2.x - p0.x
    dy2 = p2.y - p0.y
    if dx1*dy2 > dy1*dx2:
        return 1
    if dx1*dy2 < dy1*dx2:
        return -1
    if dx1*dx2 < 0 or dy1*dy2 < 0:
        return -1
    if dx1*dx1 + dy1*dy1 < dx2*dx2 + dy2*dy2:
        return 1
    return 0

def pseudo_angle(p1, p2):
    """
    @param p1: an object with x and y members
    @param p2: an object with x and y members
    @return: a number between 0 and 360 that has the same order properties as an angle
    """
    dx = p2.x - p1.x
    dy = p2.y - p1.y
    ax = math.fabs(dx)
    ay = math.fabs(dy)
    if ax + ay == 0:
        t = 0
    else:
        t = float(dy) / (ax + ay)
    if dx < 0:
        t = 2 - t
    elif dy < 0:
        t += 4
    return t * 90.0

def get_pivot_index(points):
    """
    This is the first step of the Graham scan.
    Get the index with the min y value.
    Break ties by returning the index with the max x value.
    @return: the index of the pivot point
    """
    augmented = [(p.y, -p.x, i) for (i, p) in enumerate(points)]
    return min(augmented)[2]

class GrahamScanError(Exception):
    pass

def graham_scan(points):
    """
    @return: a counter clockwise ordered list of points on the convex hull.
    """
    pivot_index = get_pivot_index(points)
    points[pivot_index], points[0] = points[0], points[pivot_index]
    pivot = points[0]
    augmented_points = [(pseudo_angle(pivot, p), p) for p in points]
    pseudo_angles, sorted_points = zip(*list(sorted(augmented_points)))
    if len(set(pseudo_angles)) != len(pseudo_angles):
        raise GrahamScanError('collinearity detected in the simple closed path')
    arr = [sorted_points[-1]] + list(sorted_points)
    M = 3
    i = 4
    while i < len(arr):
        while (ccw(arr[M], arr[M-1], arr[i]) >= 0):
            M -= 1
        M += 1
        arr[i], arr[M] = arr[M], arr[i]
        i += 1
    return arr[1:M+1]

def get_interior_quadrilateral(points):
    """
    This function can be used for interior point elimination.
    @return: four corners defining a convex quadrilateral.
    """
    pxpy = max((+p.x + p.y, p) for p in points)[1]
    pxny = max((+p.x - p.y, p) for p in points)[1]
    nxpy = max((-p.x + p.y, p) for p in points)[1]
    nxny = max((-p.x - p.y, p) for p in points)[1]
    return (pxpy, nxpy, nxny, pxny)

class TestCompGeom(unittest.TestCase):

    def setUp(self):
        npoints = 10
        random.seed(0)
        self.points = [Point((random.random(), random.random())) for i in range(npoints)]

    def test_get_pivot_index(self):
        points = (Point((0, 0)), Point((0, 1)), Point((5, 5)), Point((4, 0)), Point((3, 0)))
        pivot_index = get_pivot_index(points)
        self.assertEquals(pivot_index, 3)

    def test_pseudo_angle(self):
        pivot = self.points[0]
        pseudo_angle_index_pairs = []
        angle_index_pairs = []
        for i, p in enumerate(self.points):
            dx = p.x - pivot.x
            dy = p.y - pivot.y
            pseudo_angle_index_pairs.append((pseudo_angle(pivot, p), i))
            angle_index_pairs.append((math.atan2(dy, dx) % (2*math.pi), i))
        pseudo_angle_permutation = tuple(zip(*sorted(pseudo_angle_index_pairs))[1])
        angle_permutation = tuple(zip(*sorted(angle_index_pairs))[1])
        self.assertEquals(pseudo_angle_permutation, angle_permutation)

def get_random_point():
    x = random.random()
    y = random.random()
    x = (x - .5) * .4
    y = (y - .5) * .4
    return Point((x, y))

def main():
    # require cairo for demo but not for testing or importing
    import cairo
    # generate a bunch of random points
    npoints = 12
    random.seed(0)
    points = [get_random_point() for i in range(npoints)]
    # create the drawing context
    outfile = open('out.svg', 'w')
    size = (500, 500)
    width, height = size
    surface = cairo.SVGSurface(outfile, width, height)
    context = cairo.Context(surface)
    context.scale(width, height)
    context.translate(.5, .5)
    # paint the background
    context.set_source_rgb(.9, .9, .9)
    context.paint()
    context.set_source_rgb(0, 0, 0)
    # draw the points
    for point in points:
        radius = .01
        context.arc(point.x, point.y, radius, 0, math.pi*2)
        context.fill()
    # create the interior quadrilateral
    quad_points = get_interior_quadrilateral(points)
    # draw points inside the quadrilateral in a different color
    quad_polygon = ConvexPolygon(quad_points)
    for point in points:
        if quad_polygon.covers_point(point):
            radius = .01
            context.arc(point.x, point.y, radius, 0, math.pi*2)
            context.set_source_rgb(1, 0, 0)
            context.fill()
            context.set_source_rgb(0, 0, 0)
    # draw the edges of the quadrilateral
    context.move_to(*quad_points[0].to_tuple())
    for point in quad_points[1:]:
        context.line_to(*point.to_tuple())
    context.close_path()
    context.set_line_width(.005)
    context.stroke()
    # get the convex hull
    convex_hull = graham_scan(points)
    # draw the edges of the convex hull
    context.set_line_width(.001)
    context.set_source_rgb(0, 0, 1)
    context.move_to(*convex_hull[0].to_tuple())
    for point in convex_hull[1:]:
        context.line_to(*point.to_tuple())
    context.close_path()
    context.stroke()
    context.set_source_rgb(0, 0, 0)
    # initialize the calipers
    state = CaliperState(convex_hull)
    # get the rectangle with the minimum area
    best_rect = None
    min_area = -1
    while state.rotation < math.pi/2:
        r = list(state.gen_rectangle_vertices())
        height = math.hypot(r[1].y - r[0].y, r[1].x - r[0].x)
        width = math.hypot(r[2].y - r[1].y, r[2].x - r[1].x)
        area = width * height
        if min_area < 0 or area < min_area:
            min_area = area
            best_rect = copy.deepcopy(r)
        state = state.get_next_state()
    # draw the edges of the bounding rectangle
    #box = list(state.gen_rectangle_vertices())
    box = best_rect
    context.set_line_width(.001)
    context.set_source_rgb(0, 1, 0)
    context.move_to(*box[0].to_tuple())
    for point in box[1:]:
        context.line_to(*point.to_tuple())
    context.close_path()
    context.stroke()
    context.set_source_rgb(0, 0, 0)
    # clean up the drawing context
    surface.finish()
    outfile.close()

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    #parser.add_option('-v', '--verbose', action='store_true', dest='verbose', default=False)
    #parser.add_option('-o', '--output', dest='output_filename', metavar='FILE', help='output file')
    parser.add_option('--test', action='store_true', dest='test', default=False)
    options, args = parser.parse_args()
    if options.test:
        suite = unittest.TestLoader().loadTestsFromTestCase(TestCompGeom)
        unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        main()
