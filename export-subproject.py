"""
Pull a subproject out of the primordial soup of code.
"""

import argparse

import meta

def main(args):
    with open(args.manifest) as fin:
        module_names = [x.strip() for x in fin]
    meta.add_python_files(module_names, args.target)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--manifest', required=True,
            help='a module manifest filename')
    parser.add_argument('--target', required=True,
            help='files will be created in this existing directory')
    main(parser.parse_args())
