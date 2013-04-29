"""
Check an equation related to Markov process EM and to Sylvester equations.
"""

from StringIO import StringIO
import random
import math

import numpy as np
import scipy.linalg

import Form
import FormOut


def get_form():
    """
    @return: the body of a form
    """
    return [
            ]

def get_form_out():
    return FormOut.Report()

def slow_expm_frechet(A, E):
    n = A.shape[0]
    Q = np.zeros((2*n, 2*n))
    Q[:n, :n] = A
    Q[:n, n:] = E
    Q[n:, n:] = A
    Q_expm = scipy.linalg.expm(Q)
    return Q_expm[:n, n:]

def get_response_content(fs):

    # set up print options
    np.set_printoptions(
            linewidth=1000000,
            threshold=1000000,
            )
    out = StringIO()

    # define the number of states
    nstates = 5

    # sample a random time-reversible rate matrix
    distn = np.random.exponential(1, size=nstates)
    distn /= np.sum(distn)
    pre_Q = np.random.exponential(1, size=(nstates, nstates))
    pre_Q = pre_Q + pre_Q.T
    pre_Q = np.dot(pre_Q, np.diag(distn))
    Q = pre_Q - np.diag(np.sum(pre_Q, axis=1))
    if not np.allclose(
            np.dot(np.diag(distn), Q),
            np.dot(np.diag(distn), Q).T):
        raise Exception('does not meet detailed balance')

    # sample a random non-negative summary matrix
    C = np.random.exponential(1, size=(nstates, nstates))

    # compute matrices related to the symmetrization of the rate matrix
    D_half = np.diag(np.sqrt(distn))
    D_neg_half = np.diag(np.reciprocal(np.sqrt(distn)))
    Q_sym = np.dot(D_half, np.dot(Q, D_neg_half))
    w, U = scipy.linalg.eigh(Q_sym)

    # check that symmetrized Q is actually symmetric
    if not np.allclose(Q_sym, Q_sym.T):
        raise Exception('symmetrization fail')

    # check the reconstruction of Q from its diagonalization
    Q_sym_reconstructed = np.dot(U, np.dot(np.diag(w), U.T))
    if not np.allclose(Q_sym, Q_sym_reconstructed):
        raise Exception('incorrect eigh interpretation')

    # solve the sylvester equation AX + XB = G for X
    A = -np.diag(w)
    B = np.diag(w)
    #G = np.eye(nstates) - np.dot(
    foo = np.dot(
            np.dot(U.T, D_half),
            np.dot(C, np.dot(D_neg_half, U)))
    G = np.diag(np.diag(foo)) - foo
    X = scipy.linalg.solve_sylvester(A, B, G)
    t = 1.0
    Lam_expm = np.diag(np.exp(t*w))
    Lam_foo_expm = np.diag(t*np.diag(foo)*np.exp(t*w))
    print >> out, np.diag(foo)
    print >> out
    H = -(np.dot(Lam_expm, -X) - np.dot(-X, Lam_expm) - Lam_foo_expm)
    svd_expm_frechet = np.dot(
            np.dot(D_neg_half, U), 
            np.dot(H, np.dot(U.T, D_half)))

    print >> out, 'Q:'
    print >> out, Q
    print >> out
    print >> out, 'C:'
    print >> out, C
    print >> out
    print >> out, 'expm frechet (Q, C):'
    print >> out, slow_expm_frechet(Q, C)
    print >> out
    print >> out, 'G:'
    print >> out, G
    print >> out
    print >> out, 'X:'
    print >> out, X
    print >> out
    print >> out, 'should be G, if solve_sylvester worked:'
    print >> out, np.dot(A, X) + np.dot(X, B)
    print >> out
    print >> out, 'H:'
    print >> out, H
    print >> out
    print >> out, 'svd expm frechet:'
    print >> out, svd_expm_frechet
    print >> out
    print >> out

    # Do some debugging of the block diagonal inverse.
    H = -(
            np.dot(np.diag(w), -X) -
            np.dot(-X, np.diag(w)) -
            np.diag(np.diag(foo)))
    print >> out, 'block diagonal inverse debugging...'
    print >> out
    print >> out, 'foo:'
    print >> out, foo
    print >> out
    print >> out, 'H (should equal foo):'
    print >> out, H
    print >> out
    print >> out

    # Do some debugging of the EVD frechet derivative of expm.
    H = slow_expm_frechet(np.diag(w), foo)
    svd_expm_frechet = np.dot(
            np.dot(D_neg_half, U), 
            np.dot(H, np.dot(U.T, D_half)))
    print >> out, 'some debugging...'
    print >> out
    print >> out, 'H:'
    print >> out, H
    print >> out
    print >> out, 'svd expm frechet:'
    print >> out, svd_expm_frechet
    print >> out
    print >> out

    # Try this a different way.
    # We know the matrix J of the Jordan form of the compound matrix.
    # Try solving XJ = MX where J is this jordan part
    # and M is the compound matrix.
    n = nstates
    J = np.vstack([
        np.hstack([Lam_expm, t*np.identity(n)]),
        np.hstack([np.zeros((n,n)), Lam_expm])])
    M = np.vstack([
        np.hstack([Q, C]),
        np.hstack([np.zeros((n,n)), Q])])
    # solve the sylvester equation AX + XB = G for X
    A = -M
    B = J
    G = np.zeros((2*n,2*n))
    X = scipy.linalg.solve_sylvester(A, B, G)

    print >> out, "let's try something else..."
    print >> out
    print >> out, 'J:'
    print >> out, J
    print >> out
    print >> out, 'M:'
    print >> out, M
    print >> out
    print >> out, 'X:'
    print >> out, X
    print >> out
    print >> out, '-MX + XJ:'
    print >> out, -np.dot(M, X) + np.dot(X, J)
    print >> out

    # Solve the sylvester equation QX - XQ = I.
    # The scipy linalg solves AX + XB = Q given (a, b, q).
    #A = Q
    #B = -np.array(Q)
    #C = Q - np.diag(np.diag(Q))

    #A = np.random.exponential(1, size=(nstates, nstates))
    #B = np.random.exponential(1, size=(nstates, nstates))
    #C = np.random.exponential(1, size=(nstates, nstates))

    #print >> out, C

    #C = -np.diag(np.diag(Q))
    #I = np.eye(nstates, dtype=float)
    #X = scipy.linalg.solve_sylvester(A, B, C)

    # report stuff
    """
    print >> out, 'A:'
    print >> out, A
    print >> out
    print >> out, 'B:'
    print >> out, B
    print >> out
    print >> out, 'C:'
    print >> out, C
    print >> out
    print >> out, 'X:'
    print >> out, X
    print >> out
    print >> out, 'AX + XB:'
    print >> out, np.dot(A, X) + np.dot(X, B)
    print >> out
    """

    # show the result
    return out.getvalue()

