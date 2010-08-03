"""Convert a .hud file to a .snp file.

The input file gives a binary character vector for each OTU.
The output file is in Eigenstrat format.
Note that only the first line of the input file is used.
"""

from StringIO import StringIO
import sys

import argparse

from SnippetUtil import HandlingError
import Form
import FormOut
import Util
import hud

g_tags = ['ztools:convert']

g_default_string = """
foo 1 1 1
bar 1 1 1
baz 1 0 1
""".strip()

def process(lines):
    """
    @param lines: lines of a .hud file
    """
    names, data = hud.decode(lines)
    out = StringIO()
    for i, genotype in enumerate(data[0]):
        name = 'SNP_%d' % i
        chromosome = '1'
        morgans = '0.0'
        bases = i+1
        row = [name, chromosome, morgans, bases]
        print >> out, '\t'.join(str(x) for x in row)
    return out.getvalue().rstrip()

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('data',
                'a list of OTUs names and binary character vectors',
                g_default_string),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.EigenstratSnp('out')

def get_response_content(fs):
    return process(fs.data.splitlines()) + '\n'

def main(args):
    print process(sys.stdin)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    main(parser.parse_args())
