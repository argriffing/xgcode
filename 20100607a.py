"""Compute the Tracy-Widom statistic given a count matrix.

Following Patterson et al. each row of the input .hud file
gives an OTU name followed by presence or absence of each SNP.
"""

from StringIO import StringIO
import math
import os

import numpy as np
import argparse

from SnippetUtil import HandlingError
import Form
import Carbone


g_default_hud_string = """
IC31 1 1 0 0
IC32 1 1 1 0
IC33 1 0 1 1
IC34 0 0 1 0
""".strip()

def my_eigh(M):
    """
    Unlike numpy.eigh, the results of this function are ordered.
    The returned eigenvalues and eigenvectors are sorted
    by decreasing eigenvalue.
    The eigenvectors are returned as a list of one dimensional numpy arrays.
    If treated as a 2D array, the returned eigenvector list is the transpose
    of the eigenvector result of numpy.eig.
    The eigenvectors are normalized unit right eigenvectors.
    @param M: a symmetric numpy 2D array
    @return: eigenvalues and eigenvectors
    """
    w, vt = np.linalg.eigh(M)
    ws, vs = w.tolist(), vt.T.tolist()
    sorted_pairs = list(reversed(sorted(zip(ws, vs))))
    w, v = zip(*sorted_pairs)
    return np.array(w), [np.array(x) for x in v]

def get_tracy_widom_statistic(m, n, L):
    """
    The interpretation of the parameters is purposely vague.
    It depends on whether you are doing linkage correction.
    @param m: like the number of OTUs
    @param n: like the number of SNPs
    @param L: like a normalized principal eigenvalue
    @return: the Tracy-Widom statistic
    """
    alpha = math.sqrt(n-1) + math.sqrt(m)
    mu = (alpha*alpha) / n
    sigma = (alpha / n) * (1/math.sqrt(n-1) + 1/math.sqrt(m))**(1./3.)
    return (L - mu) / sigma

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
    evals, evecs = my_eigh(X)
    L1 = evals.sum()
    L2 = np.dot(evals, evals)
    proportion = evals[0] / L1
    # compute the relative size of the first eigenvalue
    L = m*proportion
    # compute the Tracy-Widom statistic
    x = get_tracy_widom_statistic(m, n, L)
    # do linkage correction
    n_prime = ((m+1)*L1*L1) / ((m-1)*L2 - L1*L1)
    L_prime = (m-1)*proportion
    x_prime = get_tracy_widom_statistic(m, n_prime, L_prime)
    # print some infos
    print >> out, 'number of isolates:'
    print >> out, m_full
    print >> out
    print >> out, 'total number of SNPs:'
    print >> out, n_full
    print >> out
    print >> out, 'number of informative SNPs:'
    print >> out, n
    print >> out
    print >> out, 'effective number of linkage-corrected SNPs:'
    print >> out, n_prime
    print >> out
    print >> out, 'Tracy-Widom statistic (linkage-naive):'
    print >> out, x
    print >> out
    print >> out, 'Tracy-Widom statistic (linkage-corrected):'
    print >> out, x_prime
    print >> out
    print >> out, 'proportion of variance explained by principal axis:'
    print >> out, proportion
    print >> out
    print >> out, 'eigenvalues:'
    for w in evals:
        print >> out, w
    print >> out
    print >> out, 'principal axis projection:'
    for loading, name in sorted(zip(evecs[0] * evals[0], names)):
        print >> out, '\t'.join([name, str(loading)])
    return out.getvalue().rstrip()

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('hud',
                'contents of a .hud file',
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
    parser.add_argument('--hud', required=True,
            help='a .hud file')
    main(parser.parse_args())
