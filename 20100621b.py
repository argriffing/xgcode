"""Convert a phylip alignment to a binary character alignment.

In cmdline mode the phylip is read from stdin.
"""

import sys

import argparse

from SnippetUtil import HandlingError
import Form
import Phylip
import Carbone
import hud


g_default_lines = [
        '3 20',
        'Cow       ATGGCATATCCCATACAACT',
        'Carp      ATGGCACACCCAACGCAACT',
        'Chicken   ATGGCCAACCACTCCCAACT']

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('phylip', 'non-interleaved Phylip alignment',
                '\n'.join(g_default_lines)),
            Form.ContentDisposition()]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    text = process(fs.phylip.splitlines())
    filename = 'out.hud'
    disposition = "%s; filename=%s" % (fs.contentdisposition, filename) 
    response_headers = [
            ('Content-Type', 'text/plain'),
            ('Content-Disposition', disposition)]
    return response_headers, text

def process(raw_lines):
    headers, sequences = Phylip.decode(raw_lines)
    binary_rows = Carbone.get_binary_rows(sequences)
    return hud.to_blob(headers, binary_rows)

def main(args):
    print process(sys.stdin)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    main(parser.parse_args())
