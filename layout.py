"""
Some generic layout algorithms.

For example say that I have a rectangular viewing area
which I want to populate with a bunch of panes.
In particular, some of the drawing functions
of 20110327a and 20110330a should be rewritten to use
the functions here.
"""

import unittest
import math

import numpy as np

m2pi = 2 * math.pi

class LayoutError(Exception): pass


##############################################
# This section is about angles and angle intervals.
# It currently duplicates code in SpatialTree


def get_angle(pa, pb):
    """
    @param pa: an (x, y) pair
    @param pb: an (x, y) pair
    @return: the angle towards point pb from point pa
    """
    ax, ay = pa
    bx, by = pb
    return math.atan2(by - ay, bx - ax)

class AngleInterval:
    """
    An angle interval is a (low, high) pair.
    The low and high values are each in the real number interval [0, 2*pi).
    """

    def __init__(self, low, high):
        self.low = low % m2pi
        self.high = high % m2pi

    def get_magnitude(self):
        """
        @return: the angle spanned by the interval
        """
        mag = (self.high - self.low) % m2pi
        return mag

    def get_mid_angle(self):
        """
        This is useful for drawing internal node labels.
        @return: the mid angle
        """
        mid = (self.low + self.get_magnitude() / 2.0) % m2pi
        return mid

    def update(self, other):
        """
        Modify the current interval by adding another interval.
        The union of the ranges is assumed to be contiguous
        and to span less than 2*pi radians.
        """
        triples = []
        for low in (self.low, other.low):
            for high in (self.high, other.high):
                triples.append(((high-low) % m2pi, low, high))
        best_magnitude, best_low, best_high = max(triples)
        self.low, self.high = best_low, best_high

    def contains_angle(self, theta):
        if self.low == self.high:
            return False
        phi = 0
        phi += (self.high - theta) % m2pi
        phi += (theta - self.low) % m2pi
        return phi < m2pi



##############################################
# This section is about rotating or scaling an image
# to fit inside a rectangle with a given aspect ratio.


def get_scaling_factor(current_size, max_size):
    """
    @param current_size: the width and height of the current bounding box
    @param max_size: the width and height of the target bounding box
    @return: the amount by which the current size should be scaled
    """
    cwidth, cheight = current_size
    mwidth, mheight = max_size
    xscale = float(mwidth) / float(cwidth)
    yscale = float(mheight) / float(cheight)
    return min(xscale, yscale)

def rotate_2d_origin(X, theta):
    """
    Rotation is around the origin, not around the centroid.
    @param X: each row of X is an (x, y) point
    @param theta: this is the angle of rotation
    @return: rotated points
    """
    R = np.array([
        [math.cos(theta), -math.sin(theta)],
        [math.sin(theta), math.cos(theta)]])
    return np.dot(X, R)

def rotate_2d_centroid(X, theta):
    """
    Rotation is around the centroid, not around the origin.
    @param X: each row of X is an (x, y) point
    @param theta: this is the angle of rotation
    @return: rotated points
    """
    origin = np.mean(X, axis=0)
    return rotate_2d_origin(X - origin, theta) + origin

def get_axis_aligned_size(X):
    """
    @param X: each row of X is an (x, y) point
    @return: (width, height) of bounding axis aligned rectangle
    """
    return tuple(np.max(X, axis=0) - np.min(X, axis=0))

def get_best_angle(X_in, max_size):
    """
    This is a hack.
    Search over angles in a range of 180 degrees
    to find the angle for which the scaling factor is greatest.
    @return: the best angle
    """
    sf_theta_pairs = []
    increment_count = 60
    for i in range(increment_count):
        # define the current angle
        theta = i * math.pi / increment_count
        # Get the width and height of the axis aligned bounding box
        # of the centroid-rotated points.
        X = rotate_2d_centroid(X_in, theta)
        width, height = get_axis_aligned_size(X)
        # this is a hack for zero area rectangles
        if width == 0:
            width = height / 100
        if height == 0:
            height = width / 100
        # get the scaling factor
        sf = get_scaling_factor((width, height), max_size)
        sf_theta_pairs.append((sf, theta))
    best_scale, best_angle = max(sf_theta_pairs)
    return best_angle


##############################################
# This next section is about tiling panes into a rectangular viewing area.
# The rows or columns of the tiled panes may be ragged.
# Everything in this section is about integers.


def rect_is_minimal(w, h, n):
    """
    @param w: positive integer width
    @param h: positive integer height
    @param n: positive integer minimal area
    @return: True if (w, h) satisfies constraints relating to n.
    """
    return w*h >= n and w*(h-1) < n and (w-1)*h < n

def get_rect_sizes_slow(n):
    """
    Get the minimal rectangles with area at least n.
    This is a slow O(N) algorithm.
    Each rectangle is a (w, h) pair 
    where both w and h are positive integers.
    The rectangles satisfy three constraints.
    The first constraint is that w*h >= n.
    The second constraint is that (w-1)*h < n.
    The third constraint is that w*(h-1) < n.
    @param n: a positive integer area
    @return: a collection of (width, height) pairs
    """
    rects = []
    for w in range(1, n+1):
        # For the given width,
        # find the smallest height that gives at least the target area.
        h, r = divmod(n, w)
        if r:
            h += 1
        # Append the rectangle if it is minimal.
        if rect_is_minimal(w, h, n):
            rects.append((w, h))
    return rects

def get_rect_sizes(n):
    """
    Get the minimal rectangles with area at least n.
    This is a faster O(sqrt(N)) algorithm.
    @param n: a positive integer area
    @return: a collection of (width, height) pairs
    """
    sqrt_ceil = int(math.ceil(math.sqrt(n)))
    # Get the pairs where the width is less than the sqrt ceiling.
    A = []
    for w in range(1, sqrt_ceil):
        # For the given width,
        # find the smallest height that gives at least the target area.
        h, r = divmod(n, w)
        if r:
            h += 1
        # Append the rectangle if it is minimal.
        if rect_is_minimal(w, h, n):
            A.append((w, h))
    # Get a square if it is minimal.
    if rect_is_minimal(sqrt_ceil, sqrt_ceil, n):
        B = [(sqrt_ceil, sqrt_ceil)]
    else:
        B = []
    # Get the pairs where the height is less than the sqrt ceiling.
    C = [(b, a) for a, b in A]
    # Return the list of minimal rectangles.
    return A + B + C

def min_rect_to_row_major(w, h, n):
    """
    Assist ragged rectangular layout.
    @param w: width of a minimal rectangle
    @param h: height of a minimal rectangle
    @param n: the number of items
    @return: a sequence of n (row, col) pairs
    """
    if not rect_is_minimal(w, h, n):
        raise LayoutError('expected a minimal rectangle')
    return [divmod(i, w) for i in range(n)]

def min_rect_to_col_major(w, h, n):
    """
    Assist ragged rectangular layout.
    @param w: width of a minimal rectangle
    @param h: height of a minimal rectangle
    @param n: the number of items
    @return: a sequence of (row, col) pairs
    """
    if not rect_is_minimal(w, h, n):
        raise LayoutError('expected a minimal rectangle')
    return [(b, a) for a, b in min_rect_to_row_major(h, w, n)]

def min_rect_to_row_major_matrix(w, h, n):
    """
    This is an alternative representation of the layout.
    Instead of returning a sequence of n (row, col) pairs,
    this function returns a matrix of index values.
    If a cell of the matrix is empty then its entry will be None.
    @param w: width of a minimal rectangle
    @param h: height of a minimal rectangle
    @param n: the number of items
    @return: a list of lists of index entries
    """
    M = [[None]*w for row in range(h)]
    row_col_pairs = min_rect_to_row_major(w, h, n)
    for i, (row, col) in enumerate(row_col_pairs):
        M[row][col] = i
    return M

def min_rect_to_col_major_matrix(w, h, n):
    """
    This is an alternative representation of the layout.
    Instead of returning a sequence of n (row, col) pairs,
    this function returns a matrix of index values.
    If a cell of the matrix is empty then its entry will be None.
    @param w: width of a minimal rectangle
    @param h: height of a minimal rectangle
    @param n: the number of items
    @return: a list of lists of index entries
    """
    M = [[None]*w for row in range(h)]
    row_col_pairs = min_rect_to_col_major(w, h, n)
    for i, (row, col) in enumerate(row_col_pairs):
        M[row][col] = i
    return M


class TestRaggedTilingLayout(unittest.TestCase):

    def test_get_rect_sizes_slow_9(self):
        observed = sorted(get_rect_sizes_slow(9))
        expected = [(1, 9), (2, 5), (3, 3), (5, 2), (9, 1)]
        self.assertEqual(observed, expected)

    def test_get_rect_sizes_slow_10(self):
        observed = sorted(get_rect_sizes_slow(10))
        expected = [(1, 10), (2, 5), (3, 4), (4, 3), (5, 2), (10, 1)]
        self.assertEqual(observed, expected)

    def test_get_rect_sizes_slow_31(self):
        observed = sorted(get_rect_sizes_slow(31))
        expected = [
                (1, 31), (2, 16), (3, 11), (4, 8), (5, 7), (6, 6),
                (7, 5), (8, 4), (11, 3), (16, 2), (31, 1)]
        self.assertEqual(observed, expected)

    def test_rect_size_compatibility(self):
        """
        Compare the fast implementation to the slow implementation.
        The implementations should give the same results.
        """
        for i in range(1, 100):
            observed = sorted(get_rect_sizes(i))
            expected = sorted(get_rect_sizes_slow(i))
            self.assertEqual(observed, expected)

    def test_min_rect_to_row_major_9(self):
        n = 9
        observed = min_rect_to_row_major(3, 3, n)
        expected = [
                (0,0), (0,1), (0,2),
                (1,0), (1,1), (1,2),
                (2,0), (2,1), (2,2)]
        self.assertEqual(observed, expected)

    def test_min_rect_to_row_major_10(self):
        n = 10
        observed = min_rect_to_row_major(4, 3, n)
        expected = [
                (0, 0), (0, 1), (0, 2), (0, 3),
                (1, 0), (1, 1), (1, 2), (1, 3),
                (2, 0), (2, 1)]
        self.assertEqual(observed, expected)

    def test_min_rect_to_col_major_9(self):
        n = 9
        observed = min_rect_to_col_major(3, 3, n)
        expected = [
                (0,0), (1,0), (2,0),
                (0,1), (1,1), (2,1),
                (0,2), (1,2), (2,2)]
        self.assertEqual(observed, expected)

    def test_min_rect_to_col_major_10(self):
        n = 10
        observed = min_rect_to_col_major(4, 3, n)
        expected = [
                (0, 0), (1, 0), (2, 0),
                (0, 1), (1, 1), (2, 1),
                (0, 2), (1, 2), (2, 2),
                (0, 3)]
        self.assertEqual(observed, expected)
 

class TestAngleInterval(unittest.TestCase):

    def test_degenerate_angle_interval(self):
        angle_interval = AngleInterval(0, 0)
        self.assertFalse(angle_interval.contains_angle(0))
        self.assertFalse(angle_interval.contains_angle(1.234))

    def test_small_angle_interval(self):
        angle_interval = AngleInterval(5.5, 1.1)
        self.assertTrue(angle_interval.contains_angle(0))
        self.assertFalse(angle_interval.contains_angle(2.2))
        self.assertFalse(angle_interval.contains_angle(4.4))


if __name__ == '__main__':
    unittest.main()

