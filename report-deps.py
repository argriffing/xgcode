"""
Report the const data dependencies detected among the selected modules.
"""

import argparse

import meta

def main(args):
    module_names = meta.get_module_names(
            args.manifest, args.create_all, args.create_tagged)
    module_deps, const_deps = meta.get_module_and_const_deps(module_names)
    print 'module deps:'
    for module_dep in module_deps:
        print module_dep
    print
    print 'const deps:'
    for const_dep in const_deps:
        print const_dep

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest',
            help='create xmls for snippets listed in this file')
    parser.add_argument('--create_all', action='store_true',
            help='create xmls for all snippets')
    parser.add_argument('--create_tagged',
            help='create xmls for snippets with this tag')
    main(parser.parse_args())
