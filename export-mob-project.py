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

def main(args):
    with open(args.manifest) as fin:
        module_names = [x.strip() for x in fin]
    meta.add_python_files(module_names, args.target)
    xml_target = os.path.join(mobyle_core, 'Local', 'Programs')
    path_to_auto = os.path.join(python_project, 'auto.py')
    mobyle.add_xml_files(module_names, path_to_auto, xml_target)
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
    parser.add_argument('--target', required=True,
            help='python files will be created in this existing directory')
    parser.add_argument('--mobyle_core', required=True,
            help='path to the Mobyle core directory')
    parser.add_argument('--deploy', action='store_true',
            help='deploy the xml files via mobdeploy')
    main(parser.parse_args())
