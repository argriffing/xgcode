"""
Create Mobyle-compatible xml interfaces.
"""

import os

import argparse

import mobyle


def main(args):
    with open(args.manifest) as fin:
        module_names = [x.strip() for x in fin]
    mobyle.add_xml_files(module_names, args.auto, args.target)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--manifest', required=True,
            help='a module manifest filename')
    parser.add_argument('--target', required=True,
            help='xml files will be created in this existing directory')
    parser.add_argument('--auto', required=True,
            help='path to auto.py')
    main(parser.parse_args())
