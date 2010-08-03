"""Look at perturbed Laplacians with null mass vectors.

The idea is that you can perturb a Laplacian matrix
by changing the diagonal such that the eigenvector
in the nullspace is anything you want.
Because the k-cut theorem of Fiedler applies to any
acyclic matrix, the perturbed Laplacian matrix will work
just as well as the Laplacian.
So I tried to force the nullspace vector to be proportional
to the mass vector to see if anything interesting would happen.
This does not seem to do anything useful.
Anyway this was an effort to solve a problem that I have since solved.
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

# This is the adjacency matrix of an asymmetric tree.
g_A = np.array([
    [0, 0, 0, 0, 1.0/1, 0],
    [0, 0, 0, 0, 1.0/2, 0],
    [0, 0, 0, 0, 0, 1.0/4],
    [0, 0, 0, 0, 0, 1.0/5],
    [1.0/1, 1.0/2, 0, 0, 0, 1.0/3],
    [0, 0, 1.0/4, 1.0/5, 1.0/3, 0]], dtype=float)

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
    # define the adjacency matrix
    A = g_A
    n = 6
    # define some mass distributions
    m_uniform = np.ones(n) / float(n)
    m_weighted = np.array([102, 102, 102, 102, 1, 1], dtype=float) / 410
    # make the response
    out = StringIO()
    # look at the eigendecomposition of -(1/2)HDH where D is the leaf distance matrix
    HSH = Euclid.edm_to_dccov(Euclid.g_D_b)
    W_HSH, VT_HSH = np.linalg.eigh(HSH)
    print >> out, 'W for -(1/2)HDH of the leaf distance matrix:'
    print >> out, W_HSH
    print >> out, 'VT for -(1/2)HDH of the leaf distance matrix:'
    print >> out, VT_HSH
    # look at the eigendecomposition of S given a degenerate mass distribution on the full tree
    m_degenerate = np.array([.25, .25, .25, .25, 0, 0])
    S = Euclid.edm_to_weighted_cross_product(Euclid.g_D_c, m_degenerate)
    W_S, VT_S = np.linalg.eigh(S)
    print >> out, 'W for -(1/2)(Xi)D(Xi)^T of the full distance matrix with degenerate masses:'
    print >> out, W_S
    print >> out, 'VT for -(1/2)(Xi)D(Xi)^T of the full distance matrix with degenerate masses:'
    print >> out, VT_S
    # look at the effects of various mass distributions on the MDS of the full tree
    for m in (m_uniform, m_weighted):
        # the mass distribution should sum to 1
        if not np.allclose(np.sum(m), 1):
            raise ValueError('masses should sum to 1')
        # to compute the perturbed laplacian matrix first get weighted sums
        v = np.dot(m, A)
        # now divide elementwise by the masses
        v /= m
        # subtract the adjacency matrix from the diagonal formed by elements of this vector
        Lp = np.diag(v) - A
        # now get the eigendecomposition of the pseudoinverse of the perturbed laplacian
        W_Lp_pinv, VT_Lp_pinv = np.linalg.eigh(np.linalg.pinv(Lp))
        # look at the eigendecomposition of the S matrix associated with the distance matrix of this tree
        D = Euclid.g_D_c
        S = Euclid.edm_to_weighted_cross_product(D, m)
        W_S, VT_S = np.linalg.eigh(S)
        print >> out, 'perturbed laplacian:'
        print >> out, Lp
        print >> out, 'm:', m
        print >> out, 'W for the pseudoinverse of the perturbed laplacian:'
        print >> out, W_Lp_pinv
        print >> out, 'VT for the pseudoinverse of the perturbed laplacian:'
        print >> out, VT_Lp_pinv
        print >> out, 'W for the cross product matrix:'
        print >> out, W_S
        print >> out, 'VT for the cross product matrix:'
        print >> out, VT_S
    return out.getvalue().strip()

def get_response_content(fs):
    return process() + '\n'

if __name__ == '__main__': 
    print process()
