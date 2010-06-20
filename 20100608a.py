"""Compute three principal components given a .hud file.

The output is an R frame.
The principal components are computed according to Patterson et al.
in their work on population structure.
"""

from StringIO import StringIO
import math
import os

import numpy as np
import argparse

from SnippetUtil import HandlingError
import Carbone
import EigUtil
import Form


g_default_hud_string = """
IC31 1 1 0 0
IC32 1 1 1 0
IC33 1 0 1 1
IC34 0 0 1 0
""".strip()

def process(args, raw_hud_lines):
    """
    @param args: user options from the web or cmdline
    @param hud_lines: raw lines of a .hud file
    @return: results in convenient text form
    """
    out = StringIO()
    # get the ordered names from the .hud file
    words = Carbone.get_words(raw_hud_lines)
    names = [w.name for w in words]
    # normalize the names of the isolates
    if args.clean_isolates:
        names = [Carbone.clean_isolate_element(x) for x in names]
    # create the floating point count matrix
    C_full = np.vstack([w.v for w in words])
    m_full, n_full = C_full.shape
    # remove invariant columns
    C = np.vstack([v for v in C_full.T if len(set(v))>1]).T
    # get the shape of the matrix
    m, n = C.shape
    # get the column means
    u = C.mean(axis=0)
    # get the centered and normalized counts matrix
    M = (C - u) / np.sqrt(u * (1 - u))
    # construct the sample covariance matrix
    X = np.dot(M, M.T) / n
    # get the eigendecomposition of the covariance matrix
    evals, evecs = EigUtil.eigh(X)
    # scale the eigenvectos by the eigenvalues
    pcs = [w*v for w, v in zip(evals, evecs)]
    # check for sufficient number of eigenvectors
    if len(evecs) < 3:
        raise HandlingError('the matrix has insufficient rank')
    # create the R frame
    headers = ('otu', 'pc1', 'pc2', 'pc3')
    print >> out, '\t'.join(headers)
    for i, name in enumerate(names):
        typed_row = [name, pcs[0][i], pcs[1][i], pcs[2][i]]
        if args.add_indices:
            typed_row = [i+1] + typed_row
        row = [str(x) for x in typed_row]
        print >> out, '\t'.join(row)
    return out.getvalue()

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('hud',
                'contents of a .hud file',
                g_default_hud_string),
            Form.CheckGroup('cleangroup', 'more options', [
                Form.CheckItem('add_indices',
                    'add row indices for R table compatibility', True),
                Form.CheckItem('clean_isolates',
                    'force first-column elements to be IC-prefixed', True)]),
            Form.ContentDisposition()]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    text = process(fs, fs.hud.splitlines())
    disposition = "%s; filename=%s" % (fs.contentdisposition, 'pc.table') 
    response_headers = [
            ('Content-Type', 'text/plain'),
            ('Content-Disposition', disposition)]
    return response_headers, text

def main(args):
    with open(os.path.expanduser(args.hud)) as fin_hud:
        print process(fin_hud)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--hud', required=True,
            help='a .hud file')
    parser.add_argument('--dont_add_indices',
            action='store_false', dest='add_indices',
            help='do not add R compatibile row indices')
    parser.add_argument('--dont_clean_isolates',
            action='store_false', dest='clean_isolates',
            help='do not force isolate names to be IC-prefixed')
    main(parser.parse_args())
