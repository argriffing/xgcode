"""Find a widely dispersed subset of binary strings.

Given N code words that are binary strings of length m,
find a subset of k words such that the hamming distance
between the closest two words is large.
Finding the optimal subset is too hard so we will
just retry a few times and pick the best one.
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
import iterutils
import Util
import Carbone

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
    words = [Carbone.Word(line) for line in lines]
    Carbone.validate_words(words)
    text = process(words, fs.nwords, fs.nchars)
    return [('Content-Type', 'text/plain')], text

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

def main(args):
    lines = Util.get_stripped_lines(sys.stdin)
    words = [Carbone.Word(line) for line in lines]
    Carbone.validate_words(words)
    print process(words, args.nwords, args.nchars)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--nwords', default=10, type=int,
            help='find a subset with this many OTUs')
    parser.add_argument('--nchars', default=10, type=int,
            help='find a subset with this many binary characters')
    args = parser.parse_args()
    main(args)
