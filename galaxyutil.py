"""
Things for galaxy.

These should be unit tested but are too galaxy-specific for meta.py.
"""

import unittest

g_useless_words = ('and', 'of', 'for', 'a', 'an')

def get_split_title(description, prefix_length=20):
    """
    Galaxy wants the description to be split up.
    The name should be the first part of the input description.
    The description should be the last part of the input description.
    @param description: the long description of the mini app
    @prefix_length: choose roughly this many letters of the prefix as a name
    """
    elements = description.split()
    k = 0
    for i, x in enumerate(elements):
        if k + len(x) > prefix_length:
            break
        k += len(x)
    # chop useless words from the end of the short description
    while elements[i-1].lower() in g_useless_words:
        i -= 1
    return ' '.join(elements[:i]), ' '.join(elements[i:])


class TestGalaxyUtil(unittest.TestCase):

    def test_split_title(self):
        description = 'A sudoku solver for testing presets and timeouts.'
        prefix_length = 20
        expected = ('A sudoku solver', 'for testing presets and timeouts.')
        observed = get_split_title(description, prefix_length)
        self.assertEqual(expected, observed)


if __name__ == '__main__':
    unittest.main()
