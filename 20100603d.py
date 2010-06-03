"""Convert an .ind file to a .pheno file.

The output file is in eigenstrat format.
"""

from StringIO import StringIO
import sys
import os

import argparse

from SnippetUtil import HandlingError
import Form

g_default_ind_string = """
IC1 U   Case
IC2 U   Control
IC3 U   Control
IC4 M   Ignore
""".strip()

def process(lines):
    """
    @param lines: lines of an .ind file
    @return: the single string of a .pheno file
    """
    values = []
    for line in lines:
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
                g_default_ind_string)]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    text = process(fs.ind.splitlines())
    return [('Content-Type', 'text/plain')], text

def main(args):
    with open(os.path.expanduser(args.hud)) as fin:
        print process(fin)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ind', help='.ind file')
    main(parser.parse_args())
