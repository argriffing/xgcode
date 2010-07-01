"""
List dependencies of a particularly formatted python script.
"""

import sys

import argparse

import meta

class Dep(object):

    def __init__(self):
        self.d_local_deps = {}
        self.d_external_deps = {}

    def get_immediate_local_deps(self, module_name):
        """
        @param module_name: the name of a module
        @return: a set of module names
        """
        filename = module_name + '.py'
        with open(filename) as fin:
            paragraphs = meta.get_import_paragraphs(fin)
        tiers = meta.get_import_tiers(paragraphs)

    def get_immediate_external_deps(self, module_name):
        """
        @param module_name: the name of a module
        @return: a set of module names
        """
        filename = module_name + '.py'
        pass

    def get_transitive_local_deps(self, module_name):
        """
        @param module_name: the name of a module
        @return: a set of module names
        """
        filename = module_name + '.py'
        pass

    def get_transitive_external_deps(self, module_name):
        """
        @param module_name: the name of a module
        @return: a set of module names
        """
        filename = module_name + '.py'
        pass



def process_module(module_filename):
    with open(module_filename) as fin:
        tiered_names = meta.get_tiered_names(fin)
    print module_filename
    print tiered_names
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
