"""Convert a .hud file to a .snp file.

The input file gives a binary character vector for each OTU.
The output file is in Eigenstrat format.
Note that only the first line of the input file is used.
"""

from StringIO import StringIO
import itertools
import random
import time
import sys

import argparse
import numpy as np

from SnippetUtil import HandlingError
import Form
import Util

g_default_string = """
foo 1 1 1
bar 1 1 1
baz 1 0 1
""".strip()


class WordError(Exception): pass

class Word(object):
    def __init__(self, line):
        # store the original line
        self.line = line
        # extract the name and the binary vector
        elements = line.split()
        self.name = elements[0]
        for x in elements[1:]:
            if x not in ('0', '1'):
                msg = 'expected 0 or 1 but got ' + x
                raise WordError(msg)
        self.v = np.array([int(x) for x in elements[1:]])

def validate_words(words):
    """
    Make sure the binary vectors are the same lengths.
    @param words: word objects
    """
    lengths = set(len(w.v) for w in words)
    if len(lengths) != 1:
        msg = 'each binary vector should be the same length'
        raise WordError(msg)


def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('data',
                'a list of OTUs names and binary character vectors',
                g_default_string)]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    lines = Util.get_stripped_lines(StringIO(fs.data))
    words = [Word(line) for line in lines]
    validate_words(words)
    text = process(words[0])
    return [('Content-Type', 'text/plain')], text

def process(word):
    out = StringIO()
    for i, genotype in enumerate(word.v):
        name = 'SNP_%d' % i
        chromosome = '1'
        morgans = '0.0'
        bases = i+1
        row = [name, chromosome, morgans, bases]
        print >> out, '\t'.join(str(x) for x in row)
    return out.getvalue().rstrip()

def main(args):
    lines = Util.get_stripped_lines(sys.stdin)
    words = [Word(line) for line in lines]
    validate_words(words)
    print process(words[0])

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    args = parser.parse_args()
    main(args)
