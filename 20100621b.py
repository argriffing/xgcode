"""Convert a phylip alignment to a binary character alignment.

In cmdline mode the phylip is read from stdin.
"""

import sys

import argparse

from SnippetUtil import HandlingError
import Form
import FormOut
import Phylip
import Carbone
import hud
import const

g_default_data = const.read('20100625a')

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('phylip', 'non-interleaved Phylip alignment',
                g_default_data),
            Form.RadioGroup('format_out', 'output format', [
                Form.RadioItem('hud', 'hud', True),
                Form.RadioItem('phy', 'non-interleaved phylip')]),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.Report('out.%s', ['format_out'])

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    text = process(fs, fs.phylip.splitlines())
    if fs.hud:
        filename = 'out.hud'
    elif fs.phy:
        filename = 'out.phy'
    disposition = "%s; filename=%s" % (fs.contentdisposition, filename) 
    response_headers = [
            ('Content-Type', 'text/plain'),
            ('Content-Disposition', disposition)]
    return response_headers, text

def process(fs, raw_lines):
    headers, sequences = Phylip.decode(raw_lines)
    binary_rows = Carbone.get_binary_rows(sequences)
    if fs.hud:
        return hud.encode(headers, binary_rows) + '\n'
    elif fs.phy:
        binary_seqs = [''.join(str(x) for x in row) for row in binary_rows]
        return Phylip.encode(headers, binary_seqs) + '\n'

def main(args):
    print process(args, sys.stdin)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    main(parser.parse_args())
