"""
Transform a source string into a target string using some rules. [UNFINISHED]

This uses a bidirectional graph search,
and it is meant to help look for a certain linear algebra proof.
"""

import textwrap
from StringIO import StringIO

import Form
import FormOut

g_default_rules = textwrap.dedent("""
        # r: reduce
        # b: border
        # c: center
        # i: inverse
        # p: pseudoinverse
        # s: schur complement
        # u: submatrix block

        # inverses
        ii => 1
        pp <=> 1
        i <=> iii

        # properties that are easy to prove
        cub <=> c
        cpc <=> cp
        rc <=> uc

        # properties that are harder to prove
        ibs <=> cp
        r <=> cpui
        """).strip()

g_default_initial_string = 'cps'
g_default_target_string = 'bcp'

def get_form():
    """
    @return: the body of a form
    """
    # define the form objects
    form_objects = [
            MultiLine('rules', 'rules', g_default_rules)]
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

def schur(M, nleading):
    """
    @param M: a symmetric matrix with nonsingular trailing principal submatrix
    """
    A = M[:nleading, :nleading]
    B = M[:nleading, nleading:]
    C = M[nleading:, :nleading]
    D = M[nleading:, nleading:]
    D_inv = np.linalg.inv(D)
    S = A - np.dot(B, np.dot(D_inv, C))
    return S

def schur_del(M, nkeep, ndelete):
    # compute the schur complement
    S = schur(M, nkeep)
    # return the schur complement padded with zeros
    n = nkeep + ndelete
    X = np.zeros((n, n))
    X[:nkeep, :nkeep] = S
    return X

def sample_asymmetric_matrix(n):
    B = 10.0 * np.random.rand(n, n) + 5.0
    return B

def sample_symmetric_matrix(n):
    B = sample_asymmetric_matrix(n)
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

def sample_bad_nullspace(nullity, nleading, ntrailing):
    if nullity < 2:
        msg = 'the nullity is too small for the nullspace to be bad'
        raise ValueError(msg)
    nullity_a = nullity / 2
    nullity_b = nullity - nullity_a
    A = np.random.rand()
    N = 10.0 * np.random.rand(n, nullity) - 5.0

def get_response_content(fs):
    # define some dimensions
    nleading = 5
    ntrailing = 3
    nullity = 2
    n = nleading + ntrailing
    # get a random matrix assumed to be nonsingular
    M_nonsingular = sample_asymmetric_matrix(n)
    # get a random nullspace and its associated projections
    N = 10.0 * np.random.rand(n, nullity) - 5.0
    P_to_N = ndot(N, np.linalg.inv(ndot(N.T, N)), N.T)
    P_to_N_complement = np.eye(n) - P_to_N
    # get the truncated nullspace and its associated projections
    T = N[:nleading]
    P_to_T = ndot(T, np.linalg.inv(ndot(T.T, T)), T.T)
    P_to_T_complement = np.eye(nleading) - P_to_T
    # get the singular M with assumed nonsingular principal submatrix blocks
    M_singular = ndot(P_to_N_complement, M_nonsingular, P_to_N_complement)
    # get the schur complement in M and its eigendecomposition
    S = schur(M_singular, nleading)
    S_w, S_vt = np.linalg.eigh(S)
    # Get the double sided projection of the schur complement
    # onto the orthogonal complement of the truncated nullspace.
    Sp = ndot(P_to_T_complement, S, P_to_T_complement)
    Sp_w, Sp_vt = np.linalg.eigh(Sp)
    # Make a thing that is supposed to be the same as the schur complement.
    M_singular_pinv = np.linalg.pinv(M_singular)
    mystery_pinv = ndot(
            P_to_T_complement,
            M_singular_pinv[:nleading, :nleading],
            P_to_T_complement)
    mystery = np.linalg.pinv(mystery_pinv)
    # begin the output
    np.set_printoptions(linewidth=200)
    out = StringIO()
    print >> out, 'null space (N):'
    print >> out, N
    print >> out, 'schur complement (S):'
    print >> out, S
    print >> out, 'eigenvalues of S:'
    print >> out, S_w
    print >> out, 'eigenvectors of S:'
    print >> out, S_vt
    print >> out, 'double sided projection of the schur complement'
    print >> out, 'onto the complement of the truncated nullspace of M (Sp)'
    print >> out, Sp
    print >> out, 'eigenvalues of Sp:'
    print >> out, Sp_w
    print >> out, 'eigenvectors of Sp:'
    print >> out, Sp_vt
    print >> out, 'this thing that is supposed to be the schur complement:'
    print >> out, mystery
    print >> out, 'difference from the schur complement:'
    print >> out, mystery - S
    return out.getvalue()

