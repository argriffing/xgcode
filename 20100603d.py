"""Convert an .ind file to a .pheno file.

The output file is in eigenstrat format.
"""

from StringIO import StringIO
import sys
import os

import argparse

from SnippetUtil import HandlingError
import Form
import FormOut
import iterutils

g_tags = ['pca:misc']

g_default_ind_string = """
IC1 U   Case
IC2 U   Control
IC3 U   Control
IC4 M   Ignore
""".strip()

def process(raw_lines):
    """
    @param lines: lines of an .ind file
    @return: the single string of a .pheno file
    """
    values = []
    for line in iterutils.stripped_lines(raw_lines):
        name, gender, status = line.split()
        if status == 'Control':
            v = '0'
        elif status == 'Case':
            v = '1'
        elif status == 'Ignore':
            v = '9'
        else:
            msg = 'Invalid status: ' + status
            raise Exception(msg)
        values.append(v)
    return ''.join(values)

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('ind',
                'contents of an .ind file',
                g_default_ind_string),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.EigenstratPheno('out')

def get_response_content(fs):
    return process(fs.ind.splitlines()) + '\n'

def main(args):
    with open(os.path.expanduser(args.hud)) as fin:
        print process(fin)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ind', help='.ind file')
    main(parser.parse_args())
