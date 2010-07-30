"""Look at distance properties of a node-weighted Laplacian-like matrix.

Look at distance properties of a Laplacian-like matrix
for non-uniform node weights.
This Laplacian-like matrix is the cross-product matrix S from the Abdi paper.
"""


from StringIO import StringIO
import random
import time

import numpy as np
import argparse

from SnippetUtil import HandlingError
import SnippetUtil
import Form
import FormOut
import NewickIO
import FelTree
import Euclid
import TreeSampler
import Util
import const

g_D_abdi = const.read('20100730h')

def get_form():
    """
    @return: a list of form objects
    """
    # define the list of form objects
    form_objects = []
    return form_objects

def get_form_out():
    return FormOut.Report()

def process():
    """
    @return: a multi-line string that summarizes the results
    """
    np.set_printoptions(linewidth=200)
    # set up the mass distributions
    lines = Util.get_stripped_lines(g_D_abdi.splitlines())
    D = np.array([[float(x) for x in line.split()] for line in lines])
    n = len(D)
    m_uniform = np.array([.125]*8)
    m_weighted = np.array([.1, .1, .1, .1, .1, .1, .2, .2])
    # make the response
    out = StringIO()
    for m in (m_uniform, m_weighted):
        I = np.eye(n, dtype=float)
        e = np.ones(n, dtype=float)
        E = I - np.outer(e, m)
        S = (-0.5)*np.dot(E, np.dot(D, E.T))
        s = np.diag(S)
        D_prime = np.outer(s, e) + np.outer(e, s) - 2*S
        print >> out, 'To illustrate the transformation'
        print >> out, 'of the distance matrix, we will use'
        print >> out, 'the distance matrix derived from the brain scans:'
        print >> out, D
        print >> out
        print >> out, 'The elements of the mass vector m are equal to:'
        print >> out, m
        print >> out
        print >> out, 'The centering matrix is equal to:'
        print >> out, E
        print >> out
        print >> out, 'The cross product matrix is then equal to:'
        print >> out, S
        print >> out
        print >> out, "We want to show that the cross product matrix S"
        print >> out, "will give back D computed as D = se' + e's - 2S:"
        print >> out, D_prime
        print >> out
        print >> out
    return out.getvalue().strip()

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # get the response
    result_string = process()
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, result_string

if __name__ == '__main__': 
    print process()
