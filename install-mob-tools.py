"""
This script does the second part of a two-part install.
It is meant to run in an environment that does not have lxml.
The first part of the two-part installation
uses export-mob-project.py and is run in a dependency-rich environment.
"""

import shutil
import glob
import os
import argparse
import subprocess

import mobenv

def lines_to_env_dict(raw_lines):
    d = {}
    for line in raw_lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith('#'):
            continue
        row = line.split()
        if row:
            k, v = row
            d[k] = v
    return d

def install_xml_files(stage, env_info):
    xmls = glob.glob(os.path.join(stage, 'xml-files', '*.xml'))
    xml_dir = env_info.get_xml_dir()
    for xml_src in xmls:
        xml_dst = os.path.join(xml_dir, os.path.basename(xml_src))
        shutil.copy(xml_src, xml_dst)

def install_python_subtree(stage, env_info):
    ztools_dir = os.path.dirname(env_info.auto_path)
    shutil.copytree(os.path.join(stage, 'python-files'), ztools_dir)

def main(args):
    stage = os.path.dirname(__file__)
    with open(os.path.join(stage, 'install-mob-tools.conf')) as fin:
        env_dict = lines_to_env_dict(fin)
    env_info = mobenv.create_environment_info(
            env_dict['auto_path'],
            env_dict['python_path'],
            env_dict['mob_core'])
    install_xml_files(stage, env_info)
    install_python_subtree(stage, env_info)
    if args.deploy:
        cmd = env_info.get_deploy_command()
        subprocess.call(cmd)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--deploy', action='store_true',
            help='deploy the xml files via mobdeploy')
    main(parser.parse_args())
