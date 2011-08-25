"""
Check an identity relating the pseudoinverse and two projections.

The pseudoinverse is the Moore-Penrose pseudoinverse
and the two projections are arbitrary orthogonal projections.
All matrices involved have the same rank.
Or maybe the two projections must satisfy the
equation PQP = QPQ.
I don't know if this is a general property of orthogonal projections.
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

def ndot(*args):
    M = args[0]
    for B in args[1:]:
        M = np.dot(M, B)
    return M

def get_p_del(nkeep, ndelete):
    """
    Get the deletion projection matrix.
    @param nkeep: the number of elements to keep
    @param ndelete: the number of elements to delete
    @return: an orthogonal projection matrix
    """
    n = nkeep + ndelete
    P = np.zeros((n, n))
    for i in range(nkeep):
        P[i,i] = 1.0
    return P

def get_p_centering(n):
    """
    @param n: the size of the projection matrix
    @return: an orthogonal projection matrix
    """
    e = np.ones(n)
    I = np.eye(n)
    P = np.outer(e, e) / np.inner(e, e)
    return I - P

def get_p_centering_partial(nkeep, ndelete):
    """
    Note that this returns an oblique projection matrix.
    @param nkeep: the number of elements to keep
    @param ndelete: the number of elements to delete
    @return: an oblique projection matrix
    """
    n = nkeep + ndelete
    e = np.ones(n)
    p = np.zeros(n)
    for i in range(nkeep):
        p[i] = 1.0
    P = np.outer(e, p) / np.inner(e, p)
    return np.diag(e) - P

def get_p_centering_partial_del(nkeep, ndelete):
    """
    @param nkeep: the number of elements to keep
    @param ndelete: the number of elements to delete
    @return: an orthogonal projection matrix
    """
    n = nkeep + ndelete
    e = np.ones(n)
    p = np.zeros(n)
    for i in range(nkeep):
        p[i] = 1.0
    P = np.outer(p, p) / np.inner(p, p)
    return np.diag(p) - P

def schur_del(M, nkeep, ndelete):
    # compute the schur complement
    A = M[:nkeep, :nkeep]
    B = M[:nkeep, nkeep:]
    C = M[nkeep:, nkeep:]
    C_inv = np.linalg.inv(C)
    S = A - np.dot(B, np.dot(C_inv, B.T))
    # return the schur complement padded with zeros
    n = nkeep + ndelete
    X = np.zeros((n, n))
    X[:nkeep, :nkeep] = S
    return X

def get_lhs(M, nkeep, ndelete):
    n = nkeep + ndelete
    P = get_p_centering_partial_del(nkeep, ndelete)
    return np.linalg.pinv(ndot(P, M, P))

def get_rhs(M, nkeep, ndelete):
    n = nkeep + ndelete
    P = get_p_centering_partial_del(nkeep, ndelete)
    R = get_p_centering_partial(nkeep, ndelete)
    H = get_p_centering(n)
    D = get_p_del(nkeep, ndelete)
    # define the target value in a few different ways
    target_a = schur_del(np.linalg.pinv(ndot(H, M, H)), nkeep, ndelete)
    target_b = np.linalg.pinv(ndot(P, M, P))
    target_c = np.linalg.pinv(ndot(D, R, M, R.T, D))
    # Try to find another way to get the target value using projections,
    # hopefully a way that uses projections in a way that can be shown
    # to be equivalent to taking a schur complement in a pseudoinverse.
    return np.linalg.pinv(ndot(D, R, M, R.T, D))

def get_lhs_also_old(M, nkeep, ndelete):
    """
    This is the lhs of a putative identity.
    """
    n = nkeep + ndelete
    R = get_deletion_projection(nkeep, ndelete)
    H = get_centering_projection(n)
    # build up the matrix
    X = M
    X = ndot(R, H, R, X, R, H, R)
    X = np.linalg.pinv(X)
    return X

def get_rhs_also_old(M, nkeep, ndelete):
    """
    This is the rhs of a putative identity.
    """
    n = nkeep + ndelete
    R = get_deletion_projection(nkeep, ndelete)
    H = get_centering_projection(n)
    # build up the matrix
    X = M
    X = ndot(R, H, X, H, R)
    X = np.linalg.pinv(X)
    X = ndot(R, H, X, H, R)
    return X

def get_lhs_old(M, nkeep, ndelete):
    """
    This is the lhs of a putative identity.
    """
    n = nkeep + ndelete
    R = get_deletion_projection(nkeep, ndelete)
    H = get_centering_projection(n)
    P = ndot(R, H, R)
    inside = ndot(P.T, M, P)
    return ndot(P, np.linalg.pinv(inside), P.T)

def get_rhs_old(M, nkeep, ndelete):
    """
    This is the rhs of a putative identity.
    """
    n = nkeep + ndelete
    R = get_deletion_projection(nkeep, ndelete)
    H = get_centering_projection(n)
    P = ndot(H, R, H)
    inside = ndot(P.T, M, P)
    return ndot(P, np.linalg.pinv(inside), P.T)

def sample_symmetric_matrix(n):
    B = 10.0 * np.random.rand(n, n) + 5.0
    return B + B.T

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

def get_response_content(fs):
    nkeep = 4
    ndelete = 6
    n = nkeep + ndelete
    # show the results
    np.set_printoptions(linewidth=200)
    out = StringIO()
    try:
        for i in range(1000):
            M = sample_symmetric_matrix(n)
            lhs = get_lhs(M, nkeep, ndelete)
            rhs = get_rhs(M, nkeep, ndelete)
            assert_named_equation(
                    (lhs, 'lhs'),
                    (rhs, 'rhs'))
    except Counterexample, e:
        print >> out, 'counterexample:'
        print >> out, str(e)
    else:
        print >> out, 'no counterexample was found'
    return out.getvalue()

