"""
Create Mobyle-compatible xml interfaces.
"""

import os

import argparse

import mobyle


def main(args):
    with open(args.manifest) as fin:
        module_names = [x.strip() for x in fin]
    for name in module_names:
        try:
            xml_content = mobyle.get_xml(args.auto, name)
        except mobyle.MobyleError as e:
            xml_content = None
            print e
        if xml_content:
            xml_filename = os.path.join(args.target, 'z' + name + '.xml')
            with open(xml_filename, 'w') as fout:
                fout.write(xml_content)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--manifest', required=True,
            help='a module manifest filename')
    parser.add_argument('--target', required=True,
            help='files will be created in this existing directory')
    parser.add_argument('--auto', required=True,
            help='path to auto.py')
    main(parser.parse_args())
