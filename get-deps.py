"""
List dependencies of a particularly formatted python script.
"""

import sys

import argparse

import meta

if __name__ == '__main__':
    filename = sys.argv[1]
    with open(filename) as fin:
        print meta.get_import_paragraphs(fin)
