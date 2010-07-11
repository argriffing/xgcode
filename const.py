"""
Deal with constant data.
"""

import unittest
import os

def read(name):
    """
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
