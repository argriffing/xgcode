"""
Look at Schur complements of specially structured matrices.

There is not really anything interesting here.
"""

from StringIO import StringIO

import numpy as np
import scipy
import scipy.linalg

import Form
import FormOut

class Counterexample(Exception): pass


def get_form():
    """
    @return: the body of a form
    """
    # define the form objects
    form_objects = []
    return form_objects

def get_form_out():
    return FormOut.Report()

def bott_duffin(M):
    """
    We pretend that P_L is H.
    """
    nrows, ncols = M.shape
    if nrows != ncols:
        raise ValueError('expected a square matrix')
    e = np.ones(nrows)
    I = np.eye(nrows)
    P = np.outer(e, e) / np.inner(e, e)
    H = I - P
    return np.dot(H, np.linalg.inv(np.dot(M, H) + P))

def pinvproj(M):
    """
    This could be made more efficient using double centering.
    """
    nrows, ncols = M.shape
    if nrows != ncols:
        raise ValueError('expected a square matrix')
    e = np.ones(nrows)
    I = np.eye(nrows)
    P = np.outer(e, e) / np.inner(e, e)
    H = I - P
    HMH = np.dot(H, np.dot(M, H))
    return np.linalg.inv(HMH + P) - P

def schur(M, nsmall):
    B = M[:nsmall, nsmall:]
    C = np.linalg.inv(M[nsmall:, nsmall:])
    return M[:nsmall, :nsmall] - np.dot(B, np.dot(C, B.T))

def assert_named_equation(a_name_pair, b_name_pair):
    """
    This is a helper function to check for counterexamples.
    """
    a, a_name = a_name_pair
    b, b_name = b_name_pair
    if not np.allclose(a, b):
        out = StringIO()
        print >> out, a_name + ':'
        print >> out, a
        print >> out, b_name + ':'
        print >> out, b
        raise Counterexample(out.getvalue())

def assert_pinvproj(fs, M):
    """
    Raise a Counterexample if one is found.
    """
    pinvproj_of_sub = pinvproj(M[:fs.nsmall, :fs.nsmall])
    schur_of_pinvproj = schur(pinvproj(M), fs.nsmall)
    bottduff_of_sub = bott_duffin(M[:fs.nsmall, :fs.nsmall])
    schur_of_bottduff = schur(bott_duffin(M), fs.nsmall)
    assert_named_equation(
            (pinvproj_of_sub, 'pinvproj of sub'),
            (schur_of_pinvproj, 'schur of pinvproj'))
    assert_named_equation(
            (pinvproj_of_sub, 'pinvproj of sub'),
            (bottduff_of_sub, 'bottduff of sub'))
    assert_named_equation(
            (schur_of_pinvproj, 'schur of pinvproj'),
            (schur_of_bottduff, 'schur of bottduff'))

def double_centered(M):
    """
    This is slow.
    """
    nrows, ncols = M.shape
    if nrows != ncols:
        raise ValueError('expected a square matrix')
    e = np.ones(nrows)
    I = np.eye(nrows)
    P = np.outer(e, e) / np.inner(e, e)
    H = I - P
    return np.dot(H, np.dot(M, H))

def inverse_in_H(M):
    nrows, ncols = M.shape
    if nrows != ncols:
        raise ValueError('expected a square matrix')
    e = np.ones(nrows)
    I = np.eye(nrows)
    P = np.outer(e, e) / np.inner(e, e)
    return np.linalg.inv(M + P) - P

def assert_double_centering_identities(fs, M):
    """
    Raise a Counterexample if one is found.
    """
    A = M[:fs.nsmall, :fs.nsmall]
    HMH = double_centered(M)
    HMH_A = HMH[:fs.nsmall, :fs.nsmall]
    # check an identity between two ways to compute a centered submatrix
    HAH_direct = double_centered(A)
    HAH_indirect = double_centered(HMH_A)
    assert_named_equation(
            (HAH_direct, 'double centered submatrix'),
            (HAH_indirect, 're-centered submatrix of doubly centered matrix'))
    HAH = HAH_indirect
    # This is not true:
    # check that HAH <==inverse_in_H==> H A^-1 H
    HAH_pinv = inverse_in_H(HAH)
    H_Ainv_H = double_centered(np.linalg.inv(A))
    assert_named_equation(
            (HAH_pinv, 'inverse-in-H of HAH'),
            (H_Ainv_H, 'double centered inverse of A'))
    # check a submatrix-schur-inverse commutativity in double centering
    HMH_pinv = inverse_in_H(HMH)
    schur_of_HMH_pinv = schur(HMH_pinv, fs.nsmall)
    assert_named_equation(
            (HAH_pinv, 'inverse-in-H of HAH'),
            (schur_of_HMH_pinv, 'schur of pinv of HMH'))


def sample_edm(nbig):
    """
    There is probably a nicer way to make an EDM.
    """
    nrows = nbig
    ncols = nbig-1
    X = np.random.rand(nrows, ncols)
    M = np.zeros((nbig, nbig))
    for i in range(nbig):
        for j in range(nbig):
            d = X[i] - X[j]
            M[i, j] = np.dot(d, d)
    return M

def sample_sym(n):
    M = 10 * (np.random.rand(n, n) - 0.5)
    return M + M.T

def get_response_content(fs):
    nbig = 10
    nsmall = 4
    e = np.ones(nbig)
    I = np.eye(nbig)
    Q = np.outer(e, e) / np.inner(e, e)
    H = I - Q
    M = sample_sym(nbig)
    W = double_centered(M)
    WpQ = W + Q
    # get eigendecomposition of schur complement
    W_SC = schur(W, nsmall)
    W_SC_w, W_SC_vt = scipy.linalg.eigh(W_SC)
    # get eigendecomposition of schur complement
    WpQ_SC = schur(WpQ, nsmall)
    WpQ_SC_w, WpQ_SC_vt = scipy.linalg.eigh(WpQ_SC)
    # show the results
    out = StringIO()
    np.set_printoptions(linewidth=200)
    print >> out, 'random symmetric matrix:'
    print >> out, M
    print >> out, 'after centering:'
    print >> out, W
    print >> out, 'after adding a constant eigenvector with unit eigenvalue:'
    print >> out, WpQ
    print >> out
    print >> out, 'eigendecomposition of Schur complement of centered matrix:'
    print >> out, W_SC_w
    print >> out, W_SC_vt
    print >> out
    print >> out, 'eigendecomposition of Schur complement of decentered matrix:'
    print >> out, WpQ_SC_w
    print >> out, WpQ_SC_vt
    print >> out
    return out.getvalue()

