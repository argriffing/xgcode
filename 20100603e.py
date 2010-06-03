"""Convert a .hud file and a MAT_pheno.txt file to an .ind file.

The .hud file provides the names of the OTUs.
The MAT_pheno.txt file provides the 'case-control' status.
The output file is in eigenstrat format.
"""

from StringIO import StringIO
import sys
import os

import argparse

from SnippetUtil import HandlingError
import Form
import iterutils
import Carbone


g_default_hud_string = """
IC31 1 1 1 0
IC32 1 1 1 0
IC33 1 0 1 0
IC34 1 0 1 0
""".strip()

g_default_matpheno_string = """
IC34    null
IC33    12
IC32    1
IC31    2
""".strip()


def process(hud_lines, matpheno_lines):
    """
    @param hud_lines: lines of a .hud file
    @param matpheno_lines: lines of a MAT_pheno.txt file
    @return: contents of an .ind file
    """
    # get the ordered names from the .hud file
    words = Carbone.get_words(hud_lines)
    names = [w.name for w in words]
    # get case and control status from the matpheno file
    cases = set()
    controls = set()
    for line in iterutils.stripped_lines(matpheno_lines):
        name, classification = line.split(None, 1)
        if classification == '1':
            cases.add(name)
        elif classification == '2':
            controls.add(name)
        elif classification in ('12', 'null'):
            # skip individuals classified like this
            pass
        else:
            msg = 'invalid MAT_pheno classification: ' + classification
            raise Exception(msg)
    # write the .ind file contents
    out = StringIO()
    for name in names:
        gender = 'U'
        classification = 'Ignore'
        if name in cases:
            classification = 'Case'
        elif name in controls:
            classification = 'Control'
        row = [name, gender, classification]
        print >> out, '\t'.join(row)
    return out.getvalue().rstrip()

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('hud',
                'contents of a .hud file',
                g_default_hud_string),
            Form.MultiLine('matpheno',
                'contents of a MAT_pheno.txt file',
                g_default_matpheno_string)]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    text = process(fs.hud.splitlines(), fs.matpheno.splitlines())
    return [('Content-Type', 'text/plain')], text

def main(args):
    with open(os.path.expanduser(args.hud)) as fin_hud:
        with open(os.path.expanduser(args.matpheno)) as fin_matpheno:
            print process(fin_hud, fin_matpheno)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--hud', required=True,
            help='a .hud file')
    parser.add_argument('--matpheno', required=True,
            help = 'a MAT_pheno.txt file')
    main(parser.parse_args())






