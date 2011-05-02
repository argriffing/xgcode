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

class LayoutError(Exception): pass

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
    @return: a sequence of (row, col) pairs
    """
    if not rect_is_minimal(w, h, n):
        raise LayoutError('expected a minimal rectangle')
    ngaps = w*h - n
    points = []
    for row in range(h):
        index = (h-1) - row
        if ngaps > index:
            ncols = w-1
        else:
            ncols = w
        for col in range(ncols):
            points.append((row, col))
    return points

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
    ngaps = w*h - n
    points = []
    for col in range(w):
        index = (w-1) - col
        if ngaps > index:
            nrows = h-1
        else:
            nrows = h
        for row in range(nrows):
            points.append((row, col))
    return points


class TestLayout(unittest.TestCase):

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
                (0,0), (0,1), (0,2), (0,3),
                (1,0), (1,1), (1,2),
                (2,0), (2,1), (2,2)]
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
                (0,0), (1,0), (2,0),
                (0,1), (1,1), (2,1),
                (0,2), (1,2),
                (0,3), (1,3)]
        self.assertEqual(observed, expected)


if __name__ == '__main__':
    unittest.main()

