"""
This module is for ad hoc linear clustering for visualizations.
"""

import unittest
import itertools

import Util

def gen_transitive_intervals(values, epsilon):
    """
    Create one dimensional clusters.
    If two values are closer than epsilon they will be grouped together.
    @param values: a set of numeric values
    @param epsilon: spacing for clusters
    @return: an ordered list of disjoint (low, high) inclusive intervals
    """
    if not values:
        return
    low = None
    high = None
    for value in sorted(values):
        if low is None:
            low = value
        if high is None or abs(value - high) < epsilon:
            high = value
        else:
            yield [low, high]
            low = value
            high = value
    yield [low, high]


class Discretizer:

    def __init__(self, values, max_categories, epsilon=1e-9):
        """
        If two values are closer than epsilon they will be grouped together.
        @param values: a set of numeric values
        @param max_categories: the maximum number of categories
        @param epsilon: spacing for clusters
        """
        if not values:
            raise ValueError('no values were provided')
        if max_categories < 1:
            raise ValueError('invalid max categories: %s' % max_categories)
        intervals = list(gen_transitive_intervals(values, epsilon))
        ngroups = min(max_categories, len(intervals))
        if ngroups == 1:
            self.boundaries = [[min(values), max(values)]]
        else:
            # create interval groups and value groups
            i_groups = list(Util.chopped_bresenham(intervals, ngroups))
            v_groups = [
                    list(itertools.chain.from_iterable(g)) for g in i_groups]
            self.boundaries = [(min(g), max(g)) for g in v_groups]

    def get_category_count(self):
        return len(self.boundaries)

    def get_boundaries(self):
        return self.boundaries

    def get_category_index(self, value):
        """
        @param value: a value belonging to one of the categories
        @return: the index of the category to which the value belongs
        """
        for i, (low, high) in enumerate(self.boundaries):
            if low <= value <= high:
                return i
        raise ValueError('the value %s was not found in any category' % value)


class TestDiscretizer(unittest.TestCase):

    def test_gen_transitive_intervals(self):
        # ---
        eps = .5
        values = [2, 1, 3]
        expected = [[1, 1], [2, 2], [3, 3]]
        observed = list(gen_transitive_intervals(values, eps))
        self.assertEquals(observed, expected)
        # ---
        eps = 10
        values = [3, 21, 4, 29, 18]
        expected = [[3, 4], [18, 29]]
        observed = list(gen_transitive_intervals(values, eps))
        self.assertEquals(observed, expected)

    def assert_category_count(self, values, max_categories, expected_categories):
        """
        @param values: a set of numeric values
        @param max_categories: the maximum number of categories
        @param expected_categories: the number of categories expected to be created by the discretizer
        """
        d = Discretizer(values, max_categories)
        self.assertEquals(d.get_category_count(), expected_categories)

    def test_discretizer_category_counts(self):
        self.assert_category_count([1], 1, 1)
        self.assert_category_count([1], 2, 1)
        self.assert_category_count([1, 1, 1], 1, 1)
        self.assert_category_count([1, 1, 1], 2, 1)
        self.assert_category_count([1, 1, 1], 3, 1)
        self.assert_category_count([1, 1, 1], 4, 1)
        self.assert_category_count([1, 2], 1, 1)
        self.assert_category_count([1, 2], 2, 2)
        self.assert_category_count([1, 2], 3, 2)

    def assert_boundary_integrity(self, discretizer):
        """
        @param discretizer: a L{Discretizer} object
        """
        boundaries = discretizer.get_boundaries()
        if len(boundaries) >= 2:
            for (first_low, first_high), (second_low, second_high) in zip(boundaries[:-1], boundaries[1:]):
                self.failUnless(first_high < second_low)

    def test_discretizer_boundaries(self):
        self.assert_boundary_integrity(Discretizer([1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3], 2))

    def test_discretizer_exceptions(self):
        self.assertRaises(ValueError, Discretizer, [], 3)
        self.assertRaises(ValueError, Discretizer, [1, 2, 3], 0)
        d = Discretizer([1, 2, 3], 3)
        self.assertRaises(ValueError, d.get_category_index, 4)


def main(max_categories):
    import sys
    values = [float(line) for line in sys.stdin]
    d = Discretizer(values, max_categories)
    print d.get_boundaries()

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    #parser.add_option('-v', '--verbose', action='store_true', dest='verbose', default=False)
    #parser.add_option('-o', '--output', dest='output_filename', metavar='FILE', default='-', help='output file')
    parser.add_option('--mcat', dest='max_categories', default='3', help='maximum number of categories')
    parser.add_option('--test', action='store_true', dest='test', default=False, help='run some unit tests')
    options, args = parser.parse_args()
    # run a test or run a demo
    if options.test:
        suite = unittest.TestLoader().loadTestsFromTestCase(TestDiscretizer)
        unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        max_categories = int(options.max_categories)
        main(max_categories)

