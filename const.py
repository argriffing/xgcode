"""
Deal with constant data.
"""

import unittest
import os

def read_line(name):
    """
    Return the first line stripped of leading and trailing whitespace.
    If the first line is empty it is still returned.
    @param name: something like '20480101c'
    @return: the first line of text, stripped of terminal whitespace
    """
    a, b = os.path.split(os.path.abspath(__file__))
    pathname = os.path.join(a, 'const-data', name + '.dat')
    with open(pathname) as fin:
        for line in fin:
            return line.strip()

def read(name):
    """
    Read the whole file.
    @param name: something like '20480101c'
    @return: all of the raw data
    """
    a, b = os.path.split(os.path.abspath(__file__))
    pathname = os.path.join(a, 'const-data', name + '.dat')
    with open(pathname) as fin:
        return fin.read()


class TestConst(unittest.TestCase):

    def test_path(self):
        mydir, myname = os.path.split(os.path.abspath(__file__))
        const_path = os.path.join(mydir, 'const-data')
        self.assertTrue(os.path.isdir(const_path))

if __name__ == '__main__':
    unittest.main()
