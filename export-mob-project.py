"""
Put a project onto a Mobyle server.
"""

import os
import sys
import shutil
import subprocess

import argparse

import meta
import mobyle

def get_module_names(manifest, create_all):
    """
    @param manifest: None or filename listing modules
    @param create_all: flag
    @return: module names
    """
    if bool(manifest) == bool(create_all):
        msg = 'expected exactly one of {manifest, create_all}'
        raise ValueError(msg)
    if manifest:
        with open(manifest) as fin:
            module_names = [x.strip() for x in fin]
    if create_all:
        module_names = []
        for name in os.listdir('.'):
            if re.match(r'^\d{8}[a-zA-Z]\.py$', name):
                module_name = name[:-3]
                module_names.append(module_name)
    return module_names

def main(args):
    module_names = get_module_names(args.manifest, args.create_all)
    # create the python subtree
    meta.add_python_files(module_names, args.target)
    # create the mobyle xml interface files
    xml_target = os.path.join(args.mobyle_core, 'Local', 'Programs')
    path_to_auto = os.path.join(args.target, 'auto.py')
    import_errors = mobyle.add_xml_files(
            module_names, path_to_auto, xml_target, args.short_length)
    for e in import_errors:
        print e
    # deploy the xml interface files
    path_to_mobdeploy = os.path.join(args.mobyle_core, 'Tools', 'mobdeploy')
    cmd = (path_to_mobdeploy, 'deploy')
    proc = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    proc_out, proc_err = proc.communicate()
    sys.stdout.write(proc_out)
    sys.stderr.write(proc_err)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', required=True,
            help='a module manifest filename')
    parser.add_argument('--create_all', action='store_true',
            help='create xmls for all snippets')
    parser.add_argument('--target', required=True,
            help='python files will be created in this existing directory')
    parser.add_argument('--mobyle_core', required=True,
            help='path to the Mobyle core directory')
    parser.add_argument('--deploy', action='store_true',
            help='deploy the xml files via mobdeploy')
    parser.add_argument('--short_length', type=int, default=20,
            help='max length of shortened snippet names')
    main(parser.parse_args())
