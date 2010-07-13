"""
Create Mobyle-compatible xml interfaces.
"""

import os
import re

import argparse

import mobyle


def main(args):
    if bool(args.manifest) == bool(args.create_all):
        msg = 'expected exactly one of {manifest, create_all}'
        raise ValueError(msg)
    if args.manifest:
        with open(args.manifest) as fin:
            module_names = [x.strip() for x in fin]
    if args.create_all:
        module_names = []
        for name in os.listdir('.'):
            if re.match(r'^\d{8}[a-zA-Z]\.py$', name):
                module_name = name[:-3]
                module_names.append(module_name)
    import_errors = mobyle.add_xml_files(
            module_names, args.auto, args.target, args.short_length)
    print 'import errors:'
    for e in import_errors:
        print e


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--manifest',
            help='create xmls for only snippets listed in the manifest')
    parser.add_argument('--create_all', action='store_true',
            help='create xmls for all snippets')
    parser.add_argument('--set_category',
            help='set this category for each')
    parser.add_argument('--target', required=True,
            help='xml files will be created in this existing directory')
    parser.add_argument('--auto', required=True,
            help='path to auto.py')
    parser.add_argument('--short_length', type=int, default=20,
            help='max length of shortened snippet names')
    main(parser.parse_args())
