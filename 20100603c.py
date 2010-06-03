"""Convert a .hud genotype file to a .geno genotype file.

The output file is in eigenstrat format;
it has one line per SNP
and each line has one character per individual.
"""

from StringIO import StringIO
import sys
import os

import argparse

from SnippetUtil import HandlingError
import Form
import Carbone

g_default_hud_string = """
IC1 1 1 1 0
IC2 1 1 1 0
IC3 1 0 1 0
""".strip()

def process(hud_lines):
    out = StringIO()
    words = Carbone.get_words(hud_lines)
    # for each individual get the genotype of each SNP
    array_per_individual = []
    for word in words:
        array_per_individual.append(word.v.tolist())
    # for each SNP get the genotype for each individual
    array_per_position = zip(*array_per_individual)
    for arr in array_per_position:
        print >> out, ''.join(str(x) for x in arr)
    return out.getvalue().rstrip()

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('hud',
                'a list of OTUs names and binary character vectors',
                g_default_hud_string)]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    text = process(fs.hud.splitlines())
    return [('Content-Type', 'text/plain')], text

def main(args):
    with open(os.path.expanduser(args.hud)) as fin_hud:
        print process(fin_hud)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--hud', help='.hud file')
    main(parser.parse_args())
