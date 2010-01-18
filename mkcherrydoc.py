"""
Build the documentation
and put it in CherryUtil.g_epydoc_directory.
"""

import os

import CherryUtil


if __name__ == '__main__':
    epy_cmd = ' '.join([
            'epydoc -v -o',
            CherryUtil.g_epydoc_directory,
            '*.so *.py'])
    os.system(epy_cmd)
