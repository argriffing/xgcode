"""
Pull a subproject out of the primordial soup of code.
"""

import os
import sys
import shutil

import argparse

import meta

g_default_modules = [
        'auto']

def copy_custom_modules(selected_module_names, target):
    selected_filenames = [x + '.py' for x in selected_module_names]
    for name in selected_module_names:
        source_filename = name + '.py'
        target_filename = os.path.join(target, name + '.py')
        shutil.copyfile(source_filename, target_filename)

def copy_default_modules(target):
    for name in g_default_modules:
        source_filename = name + '.py'
        target_filename = os.path.join(target, name + '.py')
        shutil.copyfile(source_filename, target_filename)

def copy_const_data(const_deps, target):
    if const_deps:
        const_data_dir = os.path.join(target, 'const-data')
        os.mkdir(const_data_dir)
        for name in const_deps:
            source_filename = os.path.join('const-data', name + '.dat')
            target_filename = os.path.join(const_data_dir, name + '.dat')
            shutil.copyfile(source_filename, target_filename)

def main(args):
    with open(args.manifest) as fin:
        module_names = [x.strip() for x in fin]
    module_deps, const_deps = meta.get_module_and_const_deps(module_names)
    selected_module_names = set(module_names) | module_deps[2]
    copy_custom_modules(selected_module_names, args.target)
    copy_default_modules(args.target)
    copy_const_data(const_deps, args.target)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--manifest',
            help='a module manifest filename')
    parser.add_argument('--target',
            help='files will be created in this existing directory')
    main(parser.parse_args())
