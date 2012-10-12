#
# Authors: Travis Oliphant, March 2002
#          Anthony Scopatz, August 2012 (Sparse Updates)
#          Jake Vanderplas, August 2012 (Sparse Updates)
#

import math

import numpy as np
from scipy import linalg


# test matrices
# http://people.sc.fsu.edu/~jburkardt/m_src/test_matrix_exponential/mexp_a.m
g_test_matrices = [
        np.array([
            [1, 0],
            [0, 2],
            ], dtype=float),
        np.array([
            [1, 3],
            [3, 2],
            ], dtype=float),
        np.array([
            [0, 1],
            [-39, -40],
            ], dtype=float),
        np.array([
            [-49, 24],
            [-64, 31],
            ], dtype=float),
        ]

# expm of test matrices
# http://people.sc.fsu.edu/~jburkardt/m_src/test_matrix_exponential/mexp_expa.m
g_test_matrices_expm = [
        np.array([
            [2.718281828459046, 0],
            [0, 7.389056098930650],
            ], dtype=float),
        np.array([
            [39.322809708033859, 46.166301438885753],
            [46.166301438885768, 54.711576854329110],
            ], dtype=float),
        np.array([
            [0, 2.718281828459046],
            [1.154822E-17, 2.718281828459046],
            ], dtype=float),
        np.array([
            [-0.735759, 0.551819],
            [-1.471518, 1.103638],
            ], dtype=float),
        ]


def expm(A):
    """
    Compute the matrix exponential using Pade approximation.

    Parameters
    ----------
    A : array or sparse matrix, shape(M,M)
        2D Array or Matrix (sparse or dense) to be exponentiated

    Returns
    -------
    expA : array, shape(M,M)
        Matrix exponential of A

    References
    ----------
    N. J. Higham,
    "The Scaling and Squaring Method for the Matrix Exponential Revisited",
    SIAM. J. Matrix Anal. & Appl. 26, 1179 (2005).

    """
    n_squarings = 0
    A_L1 = linalg.norm(A,1)
    ident = np.eye(A.shape[0], A.shape[1], dtype=A.dtype)
    if A_L1 < 1.495585217958292e-002:
        U,V = _pade3(A, ident)
    elif A_L1 < 2.539398330063230e-001:
        U,V = _pade5(A, ident)
    elif A_L1 < 9.504178996162932e-001:
        U,V = _pade7(A, ident)
    elif A_L1 < 2.097847961257068e+000:
        U,V = _pade9(A, ident)
    else:
        maxnorm = 5.371920351148152
        n_squarings = max(0, int(math.ceil(math.log(A_L1 / maxnorm, 2))))
        A /= 2**n_squarings
        U, V = _pade13(A, ident)
    R = linalg.solve(-U + V, U + V)
    for i in range(n_squarings):
        R = R.dot(R)
    return R

# implementation of Pade approximations of various degree
# using the algorithm presented in [Higham 2005]
# ident is the identity matrix
def _pade3(A, ident):
    b = (120., 60., 12., 1.)
    A2 = A.dot(A)
    U = A.dot(b[3]*A2 + b[1]*ident)
    V = b[2]*A2 + b[0]*ident
    return U,V

def _pade5(A, ident):
    b = (30240., 15120., 3360., 420., 30., 1.)
    A2 = A.dot(A)
    A4 = A2.dot(A2)
    U = A.dot(b[5]*A4 + b[3]*A2 + b[1]*ident)
    V = b[4]*A4 + b[2]*A2 + b[0]*ident
    return U,V

def _pade7(A, ident):
    b = (17297280., 8648640., 1995840., 277200., 25200., 1512., 56., 1.)
    A2 = A.dot(A)
    A4 = A2.dot(A2)
    A6 = A4.dot(A2)
    U = A.dot(b[7]*A6 + b[5]*A4 + b[3]*A2 + b[1]*ident)
    V = b[6]*A6 + b[4]*A4 + b[2]*A2 + b[0]*ident
    return U,V

def _pade9(A, ident):
    b = (17643225600., 8821612800., 2075673600., 302702400., 30270240.,
                2162160., 110880., 3960., 90., 1.)
    A2 = A.dot(A)
    A4 = A2.dot(A2)
    A6 = A4.dot(A2)
    A8 = A6.dot(A2)
    U = A.dot(b[9]*A8 + b[7]*A6 + b[5]*A4 + b[3]*A2 + b[1]*ident)
    V = b[8]*A8 + b[6]*A6 + b[4]*A4 + b[2]*A2 + b[0]*ident
    return U,V

def _pade13(A, ident):
    b = (64764752532480000., 32382376266240000., 7771770303897600.,
    1187353796428800., 129060195264000., 10559470521600., 670442572800.,
    33522128640., 1323241920., 40840800., 960960., 16380., 182., 1.)
    A2 = A.dot(A)
    A4 = A2.dot(A2)
    A6 = A4.dot(A2)
    U = A.dot(A6.dot(b[13]*A6 + b[11]*A4 + b[9]*A2) + b[7]*A6 + b[5]*A4 + b[3]*A2 + b[1]*ident)
    V = A6.dot(b[12]*A6 + b[10]*A4 + b[8]*A2) + b[6]*A6 + b[4]*A4 + b[2]*A2 + b[0]*ident
    return U,V


def main():
    for f in (linalg.expm, expm):
        for M, E in zip(g_test_matrices, g_test_matrices_expm):
            print M
            print E
            print f(M)
            print '---------------------'
        print '==================================='

if __name__ == '__main__':
    main()

