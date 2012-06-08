"""
This small module generates alignment partitions for use with BEAST.

Despite its name it is not specific to BEAST.
But I cannot imagine that it is of much general interest.
"""

import unittest

import gmpy


# this is for testing
g_start_stop_pairs = (
        # 8 of length 57
        (1, 57),
        (58, 114),
        (1 + 57*2, 57*3),
        (1 + 57*3, 57*4),
        (1 + 57*4, 57*5),
        (1 + 57*5, 57*6),
        (1 + 57*6, 57*7),
        (400, 456),
        # 4 of length 57*2
        (1, 114),
        (115, 228),
        (1 + 57*2*2, 57*2*3),
        (1 + 57*2*3, 57*2*4),
        # 2 of length 57*2*2
        (1, 228),
        (229, 456),
        # 1 of length 57*2*2*2
        (1, 456),
        )

def gen_hierarchical_slices(tile_width, start_index_in, sentinel_index_in):
    """
    @param tile_width: width of the smallest tile
    @param start_index_in: index of the first column
    @param sentinel_index_in: index of the sentinel column
    """
    ncolumns = sentinel_index_in - start_index_in
    if ncolumns < 1:
        raise ValueError('bad interval')
    if ncolumns % tile_width:
        raise ValueError('the tiles should exactly cover the interval')
    if gmpy.popcount(ncolumns / tile_width) != 1:
        raise ValueError('the number of tiles should be a power of two')
    nlevels = gmpy.scan1(ncolumns / tile_width) + 1
    for i in range(nlevels):
        width = tile_width * 2**i
        ntiles = ncolumns / width
        for j in range(ntiles):
            a = start_index_in + j*width
            b = start_index_in + (j+1)*width
            yield a, b


class TestBeastTiling(unittest.TestCase):

    def test_gen_hierarchical_slices_a(self):
        observed = tuple(
                (a, b-1) for a, b in gen_hierarchical_slices(57, 1, 457))
        expected = g_start_stop_pairs
        self.assertEqual(observed, expected)

    def test_gen_hierarchical_slices_b(self):
        observed = tuple(
                (a+1, b) for a, b in gen_hierarchical_slices(57, 0, 456))
        expected = g_start_stop_pairs
        self.assertEqual(observed, expected)


if __name__ == '__main__':
    unittest.main()

