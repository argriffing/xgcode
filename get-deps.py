"""
List dependencies of a particularly formatted python script.
"""

import sys

import argparse

import meta

def process_module(module_filename):
    with open(module_filename) as fin:
        paragraphs = meta.get_import_paragraphs(fin)
    tiers = meta.get_import_tiers(paragraphs)
    print module_filename
    print 'paragraphs:'
    print paragraphs
    print 'tiers:'
    print tiers
    print 

def main(args):
    if args.module:
        process_module(module)
    if args.manifest:
        with open(args.manifest) as fin:
            names = [x.strip() for x in fin]
        for name in names:
            filename = name + '.py'
            process_module(filename)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--module',
            help='a module filename')
    parser.add_argument('--manifest',
            help='a module manifest filename')
    main(parser.parse_args())
