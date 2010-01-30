"""
Functions to explore the Drosophila genetic reference panel.
"""

import unittest

import Util

class ChromoSkimmer:
    """
    Skim a data file.
    Each line of the file is a row of whitespace separated values,
    and the first column is the name of a chromosome.
    Remember the number of rows and the names of the chromosomes.
    All rows with the same chromosome name should be in a contiguous block.
    Yield unique chromosome names as they are encountered.
    """

    def __init__(self):
        self.name_set = set()
        self.last_name = None
        self.name_list = []
        self.linecount = 0

    def skim(self, rows):
        """
        @param rows: sequences defining rows of input values
        @raise Exception: when chromosomes are not contiguous blocks
        """
        for row in rows:
            self.linecount += 1
            name = row[0]
            if name != self.last_name:
                if name in self.name_set:
                    msg = 'each chromsome should be a contiguous block'
                    raise Exception(msg)
                yield name
                self.name_set.add(name)
                self.name_list.append(name)
                self.last_name = name

def check_chromo_monotonicity(rows):
    """
    Yield row indices.
    The first value of each input row should be the chromosome name.
    The second value of each input row should be an integer position.
    @param rows: sequences defining rows of typed input values
    @raise Exception: when chromosome positions do not monotonically increase
    """
    # scan the input file for formatting
    name_to_last_pos = {}
    for i, row in enumerate(rows):
        msg = 'the first two values of each row should be name and position'
        if len(row) < 2:
            raise Exception(msg)
        name, pos = row[0], row[1]
        if type(pos) is not int:
            raise Exception('the position should be an integer')
        last_pos = name_to_last_pos.get(name, None)
        msg = 'expected strictly increasing positions per chromosome'
        if last_pos is not None:
            if last_pos >= pos:
                raise Exception(msg)
        name_to_last_pos[name] = pos
        yield i


class TestDGRP(unittest.TestCase):

    def test_foo(self):
        pass


if __name__ == '__main__':
    unittest.main()
