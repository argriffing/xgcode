"""
Check equivalence of resistance and commute distances.

Use the GraphPCA paper for reference.
"""

from StringIO import StringIO
import itertools
import math

import numpy as np
import scipy.linalg

import Form
import FormOut
import MatrixUtil
from MatrixUtil import ndot
import combobreaker
import binarytolerance


def get_form():
    """
    @return: the body of a form
    """
    return [
            ]

def get_form_out():
    return FormOut.Report()

def sample_symmetric_pre_L(n):
    B = np.exp(np.random.randn(n, n))
    return B + B.T

def get_response_content(fs):

    # set up print options
    np.set_printoptions(
            linewidth=1000000,
            threshold=1000000,
            )
    out = StringIO()

    # define the number of vertices in the graph
    n = 6

    # sample a symmetric pre-Laplacian matrix
    pre_L = sample_symmetric_pre_L(n)
    L = np.diag(np.sum(pre_L, axis=1)) - pre_L

    # compute Gower's centered distance matrix
    G = scipy.linalg.pinv(L)

    # compute the resistance distance matrix
    R = np.zeros_like(G)
    R += np.outer(np.diag(G), np.ones(n))
    R += np.outer(np.ones(n), np.diag(G))
    R -= 2 * G

    # construct a non-expm transition matrix associated to the graph
    d = np.diag(L)
    P = np.eye(n) - np.dot(np.diag(np.reciprocal(d)), L)

    # compute first passage times according to GraphPCA
    M = np.zeros_like(P)
    for i in range(n):
        for k in range(n):
            for j in range(n):
                M[i, k] += d[j] * (G[i, j] - G[i, k] - G[k, j] + G[k, k])

    # compute commute times according to GraphPCA
    N = M + M.T

    # compute the graph volume
    volume = np.sum(d)

    print >> out, 'Laplacian matrix:'
    print >> out, L
    print >> out
    print >> out, "Gower's centered matrix:"
    print >> out, G
    print >> out
    print >> out, 'Resistance distance matrix:'
    print >> out, R
    print >> out
    print >> out, 'Non-expm transition matrix associated to the graph:'
    print >> out, P
    print >> out
    print >> out, 'Expected first-passage times:'
    print >> out, M
    print >> out
    print >> out, 'Commute times divided by graph volume:'
    print >> out, N / volume
    print >> out

    # show the result
    return out.getvalue()

