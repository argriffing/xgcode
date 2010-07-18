"""
This script does the second part of a two-part install.
It is meant to run in an environment that does not have lxml.
The first part of the two-part installation
uses export-mob-project.py and is run in a dependency-rich environment.
"""

import shutil
import glob
import os


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

def install_xml_files(stage, env_dict):
    xmls = glob.glob(os.path.join(stage, 'xml-files', '*.xml'))
    xml_dir = env_dict['xml_dir']
    for xml_src in xmls:
        xml_dst = os.path.join(xml_dir, os.path.basename(xml_src))
        shutil.copy(xml_src, xml_dst)

def install_python_subtree(stage, env_dict):
    ztools_dir = os.path.dirname(env_dict['auto_path'])
    shutil.copytree(os.path.join(stage, 'python-files'), ztools_dir)

def main():
    stage = os.path.dirname(__file__)
    with open(os.path.join(stage, 'install-mob-tools.conf')) as fin:
        env_dict = lines_to_env_dict(fin)
    install_xml_files(stage, env_dict)
    install_python_subtree(stage, env_dict)

if __name__ == '__main__':
    main()
