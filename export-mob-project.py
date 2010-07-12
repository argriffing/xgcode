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
        if not os.path.isdir(const_data_dir):
            os.mkdir(const_data_dir)
        for name in const_deps:
            source_filename = os.path.join('const-data', name + '.dat')
            target_filename = os.path.join(const_data_dir, name + '.dat')
            shutil.copyfile(source_filename, target_filename)

def add_python_files(module_names, python_project):
    module_deps, const_deps = meta.get_module_and_const_deps(module_names)
    selected_module_names = set(module_names) | module_deps[2]
    copy_custom_modules(selected_module_names, python_project)
    copy_default_modules(python_project)
    copy_const_data(const_deps, python_project)

def add_xml_files(module_names, python_project, mobyle_core):
    xml_target = os.path.join(mobyle_core, 'Local', 'Programs')
    path_to_auto = path.join(python_project, 'auto.py')
    for name in module_names:
        try:
            xml_content = mobyle.get_xml(path_to_auto, name)
        except mobyle.MobyleError as e:
            xml_content = None
            print e
        if xml_content:
            xml_filename = os.path.join(xml_target, 'z' + name + '.xml')
            with open(xml_filename, 'w') as fout:
                fout.write(xml_content)

def main(args):
    with open(args.manifest) as fin:
        module_names = [x.strip() for x in fin]
    add_python_files(module_names, args.target)
    add_xml_files(module_names, args.target, args.mobyle_core)
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
    parser.add_argument('--deploy',
            help='deploy the xml files via mobdeploy')
    main(parser.parse_args())
