"""
Put a project onto a Mobyle server.

Another usage is to create an installer
for an environment without lxml.
Example usage
python export-mob-project.py
--make_installer=ztools --universal_category=ztools
--mobyle_core=/home/mobyleuser/Mobyle96
--python_path=/home/mobyleuser/install/bin/python
--manifest=manifests/carbone.manifest --target=/home/mobyleuser/ztools
"""

import re
import os
import sys
import shutil
import subprocess
import tarfile

import argparse

import meta
import mobyle
import mobenv
import Util

def create_installer(args, cat_info, env_info, module_names):
    """
    Overwrite a staging subtree of the current directory.
    The idea it to tgz this directory and send it to the target server.
    On the target server it should be unzipped and extracted,
    and the installer script should put everything in the correct
    prespecified directories.
    """
    # define the local staging directory
    stage = args.make_installer
    # remove the old staging directory if it exists
    if os.path.isdir(stage):
        shutil.rmtree(stage)
    # create some directories
    os.makedirs(os.path.join(stage, 'python-files'))
    os.makedirs(os.path.join(stage, 'xml-files'))
    # add the appropriate python files and the const data
    meta.add_python_files(module_names, os.path.join(stage, 'python-files'))
    # create the xml files
    import_errors = mobyle.add_xml_files(
            cat_info, env_info,
            module_names, args.short_length,
            os.path.join(stage, 'xml-files'))
    for e in import_errors:
        print e
    # copy the installer script
    shutil.copy('install-mob-tools.py',
            os.path.join(stage, 'install-mob-tools.py'))
    # copy the installer script dependency
    shutil.copy('mobenv.py',
            os.path.join(stage, 'mobenv.py'))
    # create the installer configuration
    with open(os.path.join(stage, 'install-mob-tools.conf'), 'w') as fout:
        print >> fout, '#', ' '.join(sys.argv)
        print >> fout, '\t'.join(['auto_path', env_info.auto_path])
        print >> fout, '\t'.join(['python_path', env_info.python_path])
        print >> fout, '\t'.join(['mob_core', env_info.mob_core])
    # use subprocess instead of tarfile to create the tgz
    cmd = ['tar', 'czvf', stage + '.tgz', stage]
    subprocess.call(cmd)

def main(args):
    # check for flag conflicts
    if args.deploy and args.make_installer:
        msg = 'the "deploy" and "make_installer" flags are incompatible'
        raise ValueError(msg)
    # initialize the mobyle category information
    cat_info = mobyle.CategoryInfo(
            args.show_io_types, args.show_tags, args.universal_category)
    # get the module names
    module_names = meta.get_module_names(
            args.manifest, args.create_all, args.create_tagged)
    # define the environment on the target server
    auto_path = os.path.join(args.target, 'auto.py')
    env_info = mobenv.EnvironmentInfo(
            auto_path, args.python_path, args.mobyle_core)
    if args.make_installer:
        create_installer(args, cat_info, env_info, module_names)
    else:
        # create the python subtree
        meta.add_python_files(module_names, args.target)
        # create the mobyle xml interface files
        import_errors = mobyle.add_xml_files(
                cat_info, env_info,
                module_names, args.short_length, env_info.get_xml_dir())
        for e in import_errors:
            print e
    if args.deploy:
        cmd = env_info.get_deploy_command()
        subprocess.call(cmd)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest',
            help='create xmls for snippets listed in this file')
    parser.add_argument('--create_all', action='store_true',
            help='create xmls for all snippets')
    parser.add_argument('--create_tagged',
            help='create xmls for snippets with this tag')
    parser.add_argument('--target', required=True,
            help='python files will be created in this directory')
    parser.add_argument('--python_path', default='python',
            help='path to the python executable')
    parser.add_argument('--mobyle_core', required=True,
            help='path to the Mobyle core directory')
    parser.add_argument('--short_length', type=int, default=20,
            help='max length of shortened snippet names')
    parser.add_argument('--show_io_types', action='store_true',
            help='add mobyle categories for input and output types')
    parser.add_argument('--show_tags', action='store_true',
            help='add mobyle categories for module tags')
    parser.add_argument('--universal_category',
            help='add this mobyle category for all modules')
    parser.add_argument('--make_installer', metavar='DIR',
            help='create a local directory to be used for a remote install')
    parser.add_argument('--deploy', action='store_true',
            help='deploy the xml files via mobdeploy')
    main(parser.parse_args())
