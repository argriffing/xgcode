"""Convert a .hud genotype file to a .geno genotype file.

The output file is in eigenstrat format;
it has one line per SNP
and each line has one character per individual.
"""

import sys
import os

import argparse

from SnippetUtil import HandlingError
import Form
import FormOut
import hud

g_tags = ['carbone_lab']

g_hud_string = """
IC1 1 1 1 0
IC2 1 1 1 0
IC3 1 0 1 0
""".strip()

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('hud',
                'a list of OTUs names and binary character vectors',
                g_hud_string),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.EigenstratGeno('out')

def get_response_content(fs):
    return process(fs.hud.splitlines()) + '\n'

def process(raw_hud_lines):
    names, data = hud.decode(raw_hud_lines)
    columns = zip(*data)
    return '\n'.join(''.join(str(x) for x in c) for c in columns)

def main(args):
    with open(os.path.expanduser(args.hud)) as fin_hud:
        print process(fin_hud)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--hud', help='.hud file')
    main(parser.parse_args())
