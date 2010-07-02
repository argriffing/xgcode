"""
List dependencies of a particularly formatted python script.
"""

import sys

import argparse

import meta

class Dep(object):

    def __init__(self):
        self.d_deps = {}

    def get_immediate_deps(self, module_name):
        """
        @param module_name: the name of a module
        @return: three sets of module names
        """
        deps = self.d_deps.get(module_name, None)
        if deps is not None:
            return deps
        filename = module_name + '.py'
        with open(filename) as fin:
            try:
                deps = meta.get_tiered_names(fin)
            except meta.MetaError as e:
                msg = 'dependency format error in %s: %s' % (filename, e)
                raise meta.MetaError(msg)
        self.d_deps[module_name] = deps
        return deps

    def get_transitive_deps(self, module_name):
        deps = [set(), set(), set()]
        visited = set()
        unexpanded = set([module_name])
        while unexpanded:
            name = unexpanded.pop()
            visited.add(name)
            next_deps = self.get_immediate_deps(name)
            deps[0].update(next_deps[0])
            deps[1].update(next_deps[1])
            next_local_deps = next_deps[2] - visited
            deps[2].update(next_local_deps)
            unexpanded.update(next_local_deps)
        return deps


def process_module(module_filename):
    with open(module_filename) as fin:
        tiered_names = meta.get_tiered_names(fin)
    print module_filename
    print tiered_names
    print

def main(args):
    depstate = Dep()
    if args.module:
        if not args.module.endswith('.py'):
            msg = 'expected the module filename to end with .py'
            raise ValueError(msg)
        module_names = [args.module[:-3]]
    if args.manifest:
        with open(args.manifest) as fin:
            module_names = [x.strip() for x in fin]
    # get transitive module dependencies
    transitive_deps = [set(), set(), set()]
    for name in module_names:
        filename = name + '.py'
        process_module(filename)
        deps = depstate.get_transitive_deps(name)
        for a, b in zip(transitive_deps, deps):
            a.update(b)
    # get const-data dependencies for all transitive module dependencies
    const_deps = set()
    for name in set(module_names) | transitive_deps[2]:
        filename = name + '.py'
        with open(filename) as fin:
            const_deps.update(meta.get_const_deps(fin))
    print 'transitive dependencies:'
    print transitive_deps
    print 'const-data dependencies:'
    print const_deps

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--module',
            help='a module filename')
    parser.add_argument('--manifest',
            help='a module manifest filename')
    main(parser.parse_args())
