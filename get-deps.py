"""
List dependencies of a particularly formatted python script.
"""

import sys

import argparse

import meta

def main(args):
    if args.module:
        if not args.module.endswith('.py'):
            msg = 'expected the module filename to end with .py'
            raise ValueError(msg)
        module_names = [args.module[:-3]]
    if args.manifest:
        with open(args.manifest) as fin:
            module_names = [x.strip() for x in fin]
    module_deps, const_deps = meta.get_module_and_const_deps(module_names)
    print 'transitive module dependencies:'
    print module_deps
    print 'transitive read only data file dependencies:'
    print const_deps

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--module',
            help='a module filename')
    parser.add_argument('--manifest',
            help='a module manifest filename')
    main(parser.parse_args())
