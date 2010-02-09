"""Fill gaps between iterated values.
Everything here should be memory efficient.
"""

import unittest

class FillerError(Exception): pass

class FillerBase:
    def __init__(self, low, high, should_fill, err_low, err_high):
        """
        @param low: an integer or None
        @param high: an integer or None
        @param should_fill: True to fill between bounds
        @param err_low: True to raise an error on a low position
        @param err_high: True to raise an error on a high position
        """
        # do validation
        if low is not None:
            if high is not None:
                if high <= low:
                    raise ValueError('invalid low/high args')
        # get the settings
        self.low = low
        self.high = high
        self.should_fill = should_fill
        self.err_low = err_low
        self.err_high = err_high
        # initialize state variables
        self.prev = None
        self.prev_inbounds = None
    def check_monotonicity(self, pos):
        if self.prev is not None:
            if pos <= self.prev:
                raise FillerError('positions should monotonically increase')
    def check_bounds(self, pos):
        """
        @raise FillerError: if pos is oob and error checking is requested
        @return: True if the position is inbounds
        """
        if self.low is not None and pos < self.low:
            if self.err_low:
                raise FillerError('a position is less than the lower bound')
            return False
        if self.high is not None and pos > self.high:
            if self.err_high:
                raise FillerError('a position is greater than the upper bound')
            return False
        return True


class FillerCounter(FillerBase):
    def __init__(self, low, high, should_fill, err_low, err_high):
        """
        @param low: an integer or None
        @param high: an integer or None
        @param should_fill: True to fill between bounds
        @param err_low: True to raise an error on a low position
        @param err_high: True to raise an error on a high position
        """
        FillerBase.__init__(self, low, high, should_fill, err_low, err_high)
        self.npositions = 0

    def fill(self, pos):
        self.check_monotonicity(pos)
        self.prev = pos
        if not self.check_bounds(pos):
            return
        if self.should_fill:
            nfiller = 0
            if self.prev_inbounds is None:
                if self.low is not None:
                    nfiller = pos - self.low
            else:
                nfiller = (pos - self.prev_inbounds) - 1
            self.npositions += nfiller
        self.npositions += 1
        self.prev_inbounds = pos

    def finish(self):
        if self.should_fill:
            if self.prev_inbounds is not None and self.high is not None:
                nfiller = self.high - self.prev_inbounds
                self.npositions += nfiller
                self.prev_inbounds = self.high


class FillerGenerator(FillerBase):
    def __init__(self, low, high, should_fill, err_low, err_high, filler):
        """
        @param low: an integer or None
        @param high: an integer or None
        @param should_fill: True to fill between bounds
        @param err_low: True to raise an error on a low position
        @param err_high: True to raise an error on a high position
        @param filler: the value to use for filler
        """
        FillerBase.__init__(self, low, high, should_fill, err_low, err_high)
        self.filler = filler

    def fill(self, pos, value):
        self.check_monotonicity(pos)
        self.prev = pos
        if not self.check_bounds(pos):
            return
        if self.should_fill:
            nfiller = 0
            if self.prev_inbounds is None:
                if self.low is not None:
                    nfiller = pos - self.low
            else:
                nfiller = (pos - self.prev_inbounds) - 1
            for i in range(nfiller):
                yield self.filler
        yield value
        self.prev_inbounds = pos

    def finish(self):
        if self.should_fill:
            if self.prev_inbounds is not None and self.high is not None:
                nfiller = self.high - self.prev_inbounds
                for i in range(nfiller):
                    yield self.filler
                self.prev_inbounds = self.high


class TestFiller(unittest.TestCase):

    def helper(self, f_counter, f_generator, positions):
        self.fc_helper(f_counter, positions)
        return self.fg_helper(f_generator, positions)

    def fc_helper(self, f_counter, positions):
        for p in positions:
            f_counter.fill(p)
        f_counter.finish()

    def fg_helper(self, f_generator, positions):
        arr = []
        for p in positions:
            arr.extend(list(f_generator.fill(p, 1)))
        arr.extend(list(f_generator.finish()))
        return arr
    
    def test_low_high_error(self):
        T = True
        self.assertRaises(ValueError, FillerCounter, 20, 20, T, T, T)
        self.assertRaises(ValueError, FillerGenerator, 20, 20, T, T, T, 9)
        self.assertRaises(ValueError, FillerCounter, 20, 8, T, T, T)
        self.assertRaises(ValueError, FillerGenerator, 20, 8, T, T, T, 9)

    def test_montonicity(self):
        should_fill = True
        err_low = False
        err_high = False
        positions_a = [1, 16, 12, 18, 50]
        positions_b = [2, 1, 12, 18, 50]
        positions_c = [1, 2, 12, 18, 6]
        for arr in (positions_a, positions_b, positions_c):
            f_c = FillerCounter(8, 20, should_fill, err_low, err_high)
            f_g = FillerGenerator(8, 20, should_fill, err_low, err_high, 9)
            self.assertRaises(FillerError, self.fc_helper, f_c, arr)
            self.assertRaises(FillerError, self.fg_helper, f_g, arr)

    def test_no_oob_error(self):
        err_low = False
        err_high = False
        positions = [1, 12, 16, 18, 50]
        e_filled = [9, 9, 9, 9, 1, 9, 9, 9, 1, 9, 1, 9, 9]
        e_unfilled = [1, 1, 1]
        for should_fill, expected in ((True, e_filled), (False, e_unfilled)):
            f_c = FillerCounter(8, 20, should_fill, err_low, err_high)
            f_g = FillerGenerator(8, 20, should_fill, err_low, err_high, 9)
            observed = self.helper(f_c, f_g, positions)
            self.assertEqual(len(expected), f_c.npositions)
            self.assertEqual(expected, observed)

    def test_oob_error(self):
        err_low = True
        err_high = True
        positions_a = [12, 16, 18, 50]
        positions_b = [1, 12, 16, 18,]
        for should_fill in (True, False):
            for arr in (positions_a, positions_b):
                f_c = FillerCounter(8, 20, should_fill, err_low, err_high)
                f_g = FillerGenerator(8, 20, should_fill, err_low, err_high, 9)
                self.assertRaises(FillerError, self.fc_helper, f_c, arr)
                self.assertRaises(FillerError, self.fg_helper, f_g, arr)

    def test_unbounded_low(self):
        err_low = True
        err_high = True
        should_fill = True
        low = None
        high = 12
        positions = [5, 7, 9]
        expected = [1, 9, 1, 9, 1, 9, 9, 9]
        fc = FillerCounter(low, high, should_fill, err_low, err_high)
        fg = FillerGenerator(low, high, should_fill, err_low, err_high, 9)
        observed = self.helper(fc, fg, positions)
        self.assertEqual(len(expected), fc.npositions)
        self.assertEqual(expected, observed)

    def test_unbounded_high(self):
        err_low = True
        err_high = True
        should_fill = True
        low = None
        high = 12
        positions = [5, 7, 9]
        expected = [1, 9, 1, 9, 1, 9, 9, 9]
        fc = FillerCounter(low, high, should_fill, err_low, err_high)
        fg = FillerGenerator(low, high, should_fill, err_low, err_high, 9)
        observed = self.helper(fc, fg, positions)
        self.assertEqual(len(expected), fc.npositions)
        self.assertEqual(expected, observed)

    def test_skip(self):
        err_low = True
        err_high = False
        should_fill = True
        low = None
        high = 12
        positions = [5, 7, 9, 20]
        expected = [1, 9, 1, 9, 1, 9, 9, 9]
        fc = FillerCounter(low, high, should_fill, err_low, err_high)
        fg = FillerGenerator(low, high, should_fill, err_low, err_high, 9)
        observed = self.helper(fc, fg, positions)
        self.assertEqual(len(expected), fc.npositions)
        self.assertEqual(expected, observed)


if __name__ == '__main__':
    unittest.main()
