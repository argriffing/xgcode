"""Convert a .hud file to a .snp file. [UNFINISHED]

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

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('words',
                'a list of OTUs names and binary character vectors',
                g_default_string),
            Form.Integer('nwords',
                'find a subset with this many OTUs', 2),
            Form.Integer('nchars',
                'find a subset with this many binary characters', 1)]
    return form_objects


def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    lines = Util.get_stripped_lines(StringIO(fs.words))
    try:
        words = [Word(line) for line in lines]
        validate_words(words)
    except WordError, e:
        raise HandlingError(e)
    text = process(words, fs.nwords, fs.nchars)
    return [('Content-Type', 'text/plain')], text


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

def words_to_matrix(words):
    """
    @param words: validated word objects
    @return: matrix where rows are OTUs and cols are binary characters
    """
    nrows = len(words)
    ncols = len(words[0].v)
    M = np.zeros((nrows, ncols), int)
    for i, w in enumerate(words):
        for j, x in enumerate(w.v):
            M[i,j] = x
    return M

def get_separation(M, row_indices, col_indices):
    """
    @param M: a binary matrix
    @param row_indices: a set of selected row indices
    @param col_indices: a set of selected column indices
    @return: the min hamming distance among pairs of selected rows
    """
    # first extract some vectors from the matrix
    vectors = []
    for i in row_indices:
        v = np.array([M[i,j] for j in col_indices])
        vectors.append(v)
    # now get the distance
    pairs = itertools.combinations(vectors, 2)
    return min(np.dot(vb - va, vb - va) for va, vb in pairs)

def get_selections(M, nrows, ncols, nseconds):
    """
    Select a set of rows and a set of columns.
    The selections should give a large separation,
    which I am defining as the minimum hamming distance
    among all pairs of selected rows.
    @param M: a binary matrix; rows are OTUs and cols are binary characters
    @param nrows: the number of requested rows
    @param ncols: the number of requested columns
    @param nseconds: a time limit for the search
    @return: a row index set and column index set
    """
    nrows_total, ncols_total = M.shape
    best_selections = None
    best_separation = None
    t = time.time()
    while time.time() - t < nseconds:
        row_indices = set(random.sample(range(nrows_total), nrows))
        col_indices = set(random.sample(range(ncols_total), ncols))
        d = get_separation(M, row_indices, col_indices)
        if (best_separation is None) or (best_separation < d):
            best_separation = d
            best_selections = (row_indices, col_indices)
    return best_selections

def process(words, nwords, nchars, nseconds=2.0):
    """
    @param words: a sequence of word objects
    @param nwords: find a subset of this many word objects
    @param nchars: find a subset of this many binary characters
    @param nseconds: a time limit
    @return: multiline string
    """
    out = StringIO()
    if len(words) < nwords:
        msg = 'the number of OTUs is smaller than the desired sample'
        raise HandlingError(msg)
    if len(words[0].v) < nchars:
        msg = 'the number of characters is smaller than the desired sample'
        raise HandlingError(msg)
    # create the matrix
    M = words_to_matrix(words)
    # select row and column indices
    row_indices, col_indices = get_selections(M, nwords, nchars, nseconds)
    sorted_row_indices = list(sorted(row_indices))
    sorted_col_indices = list(sorted(col_indices))
    # print the separation
    d = get_separation(M, row_indices, col_indices)
    print >> out, 'best separation:', d
    # print the index selections
    print >> out, 'selected row indices:', sorted_row_indices
    print >> out, 'selected column indices:', sorted_col_indices
    # print some selected values
    for i in sorted_row_indices:
        w = words[i]
        s = ' '.join(str(M[i,j]) for j in sorted_col_indices)
        print >> out, w.name + '\t' + s
    return out.getvalue().rstrip()

def process(raw_lines):
    out = StringIO()
    line = Util.get_first(Util.stripped_lines(raw_lines))
    otu_name, genotype_string = line.split(None, 1)
    genotypes = genotype_string.split()
    for i, genotype in enumerate(genotypes):
        name = 'SNP_' + str(i)
        chromosome = '1'
        morgans = '0.0'
        bases = i+1
        row = [name, chromosome, morgans, bases]
        print >> out, '\t'.join(str(x) for x in row)
    return out.getvalue().rstrip()

def main(args):
    print process(sys.stdin)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    args = parser.parse_args()
    main(args)
