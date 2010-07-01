"""
Deal with constant data.
"""

import os

def read(name):
    """
    @param name: something like '20480101c'
    @return: all of the raw data
    """
    pathname = os.path.join('const-data', name + '.dat')
    with open(pathname) as fin:
        return fin.read()
