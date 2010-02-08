"""
This module is for layouts that use monospaced fonts.
"""

import unittest


def _justify_helper(value, width, padding_value, left):
    """
    @param value: this value will be justified
    @param width: the width of the justified line
    @param padding_value: the value used as padding
    @param left: True for left justification
    @return: a justified line
    """
    padding_width = width - len(value)
    if padding_width < 0:
        raise ValueError('not enough width has been allocated for this layout')
    padding = padding_value * padding_width
    if left:
        return value + padding
    else:
        return padding + value

def left_justify(value, width, padding_value):
    """
    @param value: this value will be justified
    @param width: the width of the justified line
    @param padding_value: the value used as padding
    @return: a justified line
    """
    return _justify_helper(value, width, padding_value, True)

def right_justify(value, width, padding_value):
    """
    @param value: this value will be justified
    @param width: the width of the justified line
    @param padding_value: the value used as padding
    @return: a justified line
    """
    return _justify_helper(value, width, padding_value, False)

def _ruler_helper(label_offset_pairs, max_width, spacing_value):
    """
    @param label_offset_pairs: a list of (label_string, offset) pairs
    @param max_width: the maximum length of the string
    @param spacing_value: the string used to space the labels
    @return: a string consisting of spaced labels
    """
    arr = []
    current_width = 0
    for i, (label, offset) in enumerate(label_offset_pairs):
        # If writing this label would make the string exceed the maximum allowed width
        # then we cannot write this label.
        if offset + len(label) > max_width:
            continue
        # calculate the amount of spacing we need
        spacing_width = offset - current_width
        # If there is not sufficient spacing
        # then we cannot write this label.
        if i > 0 and spacing_width < 1:
            continue
        # add the label to the string
        s = (spacing_value * spacing_width) + label
        arr.append(s)
        current_width += len(s)
    return ''.join(arr)

def get_ruler_line(low, high):
    """
    Get a string that is a line of numbers.
    @param low: the integer lower bound of the range
    @param high: the integer upper bound of the range
    @return: a string that marks the low number and intervals of ten where possible
    """
    if low > high:
        raise ValueError('the lower bound should not exceed the upper bound')
    max_width = 1 + high - low
    label_offset_pairs = []
    for offset in range(max_width):
        tick = offset + low
        if tick == low or tick % 10 == 0:
            label_offset_pairs.append((str(tick), offset))
    spacer = ' '
    return _ruler_helper(label_offset_pairs, max_width, spacer)

def get_codon_ruler_line(low, high):
    """
    Get a string that is a line of numbers.
    @param low: the lower codon index
    @param high: the upper codon index
    @return: a string that marks the low index and high index in intervals of five where possible
    """
    if low > high:
        raise ValueError('the lower bound should not exceed the upper bound')
    ncodons = 1 + high - low
    max_width = ncodons * 3
    label_offset_pairs = []
    for codon_offset in range(ncodons):
        offset = codon_offset * 3
        tick = codon_offset + low
        if tick == low or tick % 5 == 0:
            label_offset_pairs.append((str(tick), offset))
    spacer = ' '
    return _ruler_helper(label_offset_pairs, max_width, spacer)


class TestMonospace(unittest.TestCase):

    def assert_ruler(self, low, high, expected_line):
        """
        @param low: the integer lower bound of the range
        @param high: the integer upper bound of the range
        @param expected_line: the expected line
        """
        self.assertEquals(get_ruler_line(low, high), expected_line)

    def test_ruler(self):
        self.assert_ruler(1, 23, '1        10        20')
        self.assert_ruler(1000, 1021, '1000      1010')

    def assert_codon_ruler(self, low, high, expected_line):
        """
        @param low: the lower codon index
        @param high: the upper codon index
        @param expected_line: the expected line
        """
        self.assertEquals(get_codon_ruler_line(low, high), expected_line)

    def test_codon_ruler(self):
        expected = '1           5              10             15             20'
        self.assert_codon_ruler(1, 23, expected)
        expected = '1000           1005           1010           1015           1020'
        self.assert_codon_ruler(1000, 1021, expected)

    def test_justify_string(self):
        observed = left_justify('hello world', 20, '.')
        expected = 'hello world.........'
        self.assertEquals(observed, expected)
        observed = right_justify('hello world', 20, '.')
        expected = '.........hello world'
        self.assertEquals(observed, expected)
        observed = right_justify('hello world', 11, '.')
        expected = 'hello world'
        self.assertEquals(observed, expected)

    def test_justify_list(self):
        observed = left_justify(list('hello'), 10, ['&nbsp;'])
        expected = ['h', 'e', 'l', 'l', 'o', '&nbsp;', '&nbsp;', '&nbsp;', '&nbsp;', '&nbsp;']
        self.assertEquals(observed, expected)
        observed = right_justify(list('hello'), 10, ['&nbsp;'])
        expected = ['&nbsp;', '&nbsp;', '&nbsp;', '&nbsp;', '&nbsp;', 'h', 'e', 'l', 'l', 'o']
        self.assertEquals(observed, expected)
        observed = right_justify(list('hello'), 5, ['&nbsp;'])
        expected = ['h', 'e', 'l', 'l', 'o']
        self.assertEquals(observed, expected)


def main():
    pass

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    #parser.add_option('-v', '--verbose', action='store_true', dest='verbose', default=False)
    #parser.add_option('-o', '--output', dest='output_filename', metavar='FILE', default='-', help='output file')
    parser.add_option('--test', action='store_true', dest='test', default=False, help='run some unit tests')
    options, args = parser.parse_args()
    # run a test or run a demo
    if options.test:
        suite = unittest.TestLoader().loadTestsFromTestCase(TestMonospace)
        unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        main()

