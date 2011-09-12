"""
Check a formula for the pseudoinverse of a double centered matrix.

In particular we are interested in the case where the original
matrix has zero mean of elements but is not double centered.
"""

from StringIO import StringIO

import numpy as np

import Form
import FormOut

def get_form():
    """
    @return: the body of a form
    """
    # define the form objects
    form_objects = []
    return form_objects

def get_form_out():
    return FormOut.Report()

def ndot(*args):
    M = args[0]
    for B in args[1:]:
        M = np.dot(M, B)
    return M

def double_centered(M):
    n = len(M)
    e = np.ones(n)
    I = np.eye(n)
    P = np.outer(e, e) / np.inner(e, e)
    H = I - P
    return ndot(H, M, H)

def mean_removed(M):
    return M - np.mean(M)

def sample_asymmetric_matrix(n):
    B = 10.0 * np.random.rand(n, n) - 5.0
    return B

def sample_symmetric_matrix(n):
    B = sample_asymmetric_matrix(n)
    return B + B.T

def get_response_content(fs):
    n = 5
    e = np.ones(n)
    J = np.outer(e, e)
    M = sample_asymmetric_matrix(n)
    M_sum = np.sum(M)
    MJM = ndot(M, J, M)
    R = mean_removed(M)
    RJR = ndot(R, J, R)
    # begin the output
    np.set_printoptions(linewidth=200)
    out = StringIO()
    print >> out, 'original matrix (M):'
    print >> out, M
    print >> out, 'MJM:'
    print >> out, MJM
    print >> out, 'mean-removed matrix (R):'
    print >> out, R
    print >> out, 'grand mean of R:'
    print >> out, np.mean(R)
    print >> out, 'RJR:'
    print >> out, RJR
    return out.getvalue()

