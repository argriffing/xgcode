"""
This script tries to detect whether a LaTeX package is installed.

The user specifies the name of a LaTeX package on the command line,
and this script tries to detect whether or not
the package has been installed.
This script itself is not so useful because the command line program
kpsewhich could be used instead;
the value of this python script is as a driver to demonstrate the
latexutil programmatic wrapper on top of kpsewhich.
"""

import argparse

import latexutil

def query_packages(package_names):
    names = latexutil.check_packages(package_names)
    if names:
        print 'found the following packages as .sty files:'
        for name in names:
            print name
    else:
        print 'none of these packages were found as .sty files'

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('packages', metavar='PKG', nargs='+',
            help='a LaTeX package name')
    args = parser.parse_args()
    query_packages(args.packages)
