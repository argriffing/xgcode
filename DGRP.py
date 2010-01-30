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

    def skim(self, fin):
        """
        @param fin: a file open for reading
        """
        for line in Util.stripped_lines(fin):
            self.linecount += 1
            name = line.split()[0]
            if name != self.last_name:
                if name in self.name_set:
                    msg = 'each chromsome should be a contiguous block'
                    raise Exception(msg)
                yield name
                self.name_set.add(name)
                self.name_list.append(name)
                self.last_name = name



class TestDGRP(unittest.TestCase):

    def test_foo(self):
        pass


if __name__ == '__main__':
    unittest.main()
