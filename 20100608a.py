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
from Form import RadioItem




# TODO add raw R output disposition for downloading

g_default_hud_string = """
IC31 1 1 0 0
IC32 1 1 1 0
IC33 1 0 1 1
IC34 0 0 1 0
""".strip()

def process(hud_lines):
    """
    @param hud_lines: lines of a .hud file
    @return: results in convenient text form
    """
    out = StringIO()
    # get the ordered names from the .hud file
    words = Carbone.get_words(hud_lines)
    names = [w.name for w in words]
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
    headers = ('name', 'pc1', 'pc2', 'pc3')
    print >> out, '\t'.join(headers)
    for i, name in enumerate(names):
        typed_row = (i+1, name, pcs[0][i], pcs[1][i], pcs[2][i])
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
            Form.RadioGroup('contentdisposition', 'delivery options', [ 
                RadioItem('inline', 'view', True), 
                RadioItem('attachment', 'download')])] 
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    text = process(fs.hud.splitlines())
    response_headers = [('Content-Type', 'text/plain')] 
    disposition = "%s; filename=%s" % (fs.contentdisposition, 'pc.table') 
    response_headers.append(('Content-Disposition', disposition)) 
    return response_headers, text

def main(args):
    with open(os.path.expanduser(args.hud)) as fin_hud:
        print process(fin_hud)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--hud', required=True,
            help='a .hud file')
    main(parser.parse_args())
