"""
Try to implement Algorithm 6.4 of a paper by Awad Al-Mohy and Nick Higham.

Computing the Frechet derivative of the matrix exponential,
with an application to condition number estimation.
The magic constants in this module are tuned
for double precision floating point computations.
"""

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



# Coefficients of degree 13 Pade approximant.
# From Algorithm 2.3 of Higham 2005,
# The Scaling and Squaring Method for the Matrix Exponential Revisited.
# b0 through b13
btable = np.array([
        64764752532480000, 32382376266240000, 7771770303897600,
        1187353796428800, 129060195264000, 10559470521600,
        670442572800, 33522128640, 1323241920,
        40840800, 960960, 16380, 182, 1,
        ], dtype=float)
    
b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13 = btable

table61_theta = [
        None,
        3.65e-8, 5.32e-4, 1.50e-2, 8.54e-2, 2.54e-1,
        5.41e-1, 9.50e-1, 1.47e0, 2.10e0, 2.81e0,
        3.60e0, 4.46e0, 5.37e0, 6.33e0, 7.34e0,
        8.37e0, 9.44e0, 1.05e1, 1.17e1, 1.28e1,
        ]

table61_l = [
        None,
        2.11e-8, 3.56e-4, 1.08e-2, 6.49e-2, 2.00e-1,
        4.37e-1, 7.83e-1, 1.23e0, 1.78e0, 2.42e0,
        3.13e0, 3.90e0, 4.74e0, 5.63e0, 6.56e0,
        7.52e0, 8.53e0, 9.56e0, 1.06e1, 1.17e1,
        ]


def alg34(A, E, k):
    """
    This is a naive transcription of (3.4).
    @return: M_k
    """
    if k < 1:
        raise Exception
    elif k == 1:
        return E
    else:
        return np.dot(alg34(k-1), A) + np.dot(np.linalg.matrix_power(A, k-1), E)

def alg612(A, E, m):
    """
    This is a naive transcription of (6.12) and (6.13).
    @return: L_u_m, L_v_m
    """
    Lua = np.zeros_like(A)
    Lub = np.zeros_like(A)
    Lv = np.zeros_like(A)
    for k in range(1, (m-1)//2 + 1):
        M2k = alg34(A, E, 2*k)
        Lua += btable[2*k + 1] * M2k
        Lv += btable[2*k] * M2k
    for k in range(0, (m-1)//2 + 1):
        Lub += btable[2*k + 1] * np.linalg.matrix_power(A, 2*k)
    Lu = np.dot(A, Lua) + np.dot(E, Lub)
    return Lu, Lv

def alg64_26_27(U, V, Lu, Lv):
    """
    This is a naive transcription of lines 26 and 27 of Algorithm 6.4.
    The paper has a comment that the linear systems at lines 26 and 27
    have the same coefficient matrix, so an LU factorization
    can be computed once and reused.
    """
    R = linalg.solve(-U + V, U + V)
    L = linalg.solve(-U + V, Lu + Lv + np.dot(Lu - Lv, R))
    return R, L

def alg611(X, m):
    """
    This is a naive transcription of definition (6.11).
    It is a decomposition of a sum into odd and even parts.
    """
    U = np.zeros_like(X)
    V = np.zeros_like(X)
    for k in range(0, (m-1)//2 + 1):
        X2k = np.linalg.matrix_power(X, 2*k)
        U += btable[2*k + 1] * X2k
        V += btable[2*k] * X2k
    return np.dot(X, U), V

def alg614(U, V):
    return linalg.solve(-U + V, U + V)

# from scipy
def _pade13(A, ident):
    """
    b = (64764752532480000., 32382376266240000., 7771770303897600.,
    1187353796428800., 129060195264000., 10559470521600., 670442572800.,
    33522128640., 1323241920., 40840800., 960960., 16380., 182., 1.)
    """
    b = btable
    A2 = A.dot(A)
    A4 = A2.dot(A2)
    A6 = A4.dot(A2)
    U = A.dot(A6.dot(b[13]*A6 + b[11]*A4 + b[9]*A2) + b[7]*A6 + b[5]*A4 + b[3]*A2 + b[1]*ident)
    V = A6.dot(b[12]*A6 + b[10]*A4 + b[8]*A2) + b[6]*A6 + b[4]*A4 + b[2]*A2 + b[0]*ident
    return U,V

def hardcoded_13_13_pade_approximant(A, I):
    """
    Lines 17, 18, 19, 20 of algorithm 2.3 of Higham 2005.
    """
    A2 = np.linalg.matrix_power(A, 2)
    A4 = np.linalg.matrix_power(A2, 2)
    A6 = np.dot(A4, A2)
    U = np.dot(
            A,
            np.dot(
                A6,
                b13*A6 + b11*A4 + b9*A2,
                ) + b7*A6 + b5*A4 + b3*A2 + b1*I)
    V = np.dot(
            A6,
            b12*A6 + b10*A4 + b8*A2,
            ) + b6*A6 + b4*A4 + b2*A2 + b0*I
    return U, V

def alg61(A):
    """
    Scaling and squaring algorithm for exponential.
    """
    # XXX for now avoid stopping early
    """
    for m in (3, 5, 7, 9):
        if linalg.norm(A, 1) < table61_theta[m]:
            print 'using m =', m
            U, V = alg611(A, m)
            return alg614(U, V)
    """
    m = 13
    s = int(math.ceil(math.log(linalg.norm(A, 1), 2) / table61_theta[m]))
    A /= 2**s
    """
    print 'using m =', m
    print 'using s =', s
    U, V = alg611(A, m)
    X = alg614(U, V)
    """
    #U, V = hardcoded_13_13_pade_approximant(A, np.eye(A.shape[0]))
    U, V = _pade13(A, np.eye(A.shape[0]))
    r13 = linalg.solve(-U + V, U + V)
    return np.linalg.matrix_power(r13, 2**s)
    """
    for i in range(s):
        X = np.linalg.matrix_power(X, 2)
    return X
    """

def alg64(A, E):
    """
    Scaling and squaring algorithm for exponential and Frechet derivative.
    Given A, E, this algorithm computes R and L by scaling and squaring.
    The matrices A and E are square nxn with complex entries.
    Matrix R is expm(A), and matrix L is Lexp(A, E) in their notation.
    @return: R, L
    """
    if A.ndim != 2:
        raise Exception
    n = A.shape[0]
    if A.shape != (n, n):
        raise Exception
    if E.shape != (n, n):
        raise Exception
    for m in (3, 5, 7, 9):
        if linalg.norm(A, 1) < table61_l[m]:
            # Evaluate U = um(A) and V = vm(A), using (6.11).
            U, V = alg611(A, m)
            # Evaluate Lu = Lum(A, E) and Lv = Lvm(A, E),
            # using (6.12) and (6.13).
            Lu, Lv = alg612(A, E, m)
            return alg64_26_27(U, V, Lu, Lv)
    s = int(math.ceil(math.log(linalg.norm(A, 1) / table61_l[13], 2)))
    A /= 2**s
    E /= 2**s
    A2 = np.linalg.matrix_power(A, 2)
    A4 = np.linalg.matrix_power(A2, 2)
    A6 = np.dot(A2, A4)
    M2 = np.dot(A, E) + np.dot(E, A)
    M4 = np.dot(A2, M2) + np.dot(M2, A2)
    M6 = np.dot(A4, M2) + np.dot(M4, A2)
    W1 = b13*A6 + b11*A4 + b9*A2
    I = np.eye(n)
    W2 = b7*A6 + b5*A4 + b3*A2 + b1*I
    Z1 = b12*A1 + b10*A4 + b8*A2
    Z2 = b6*A6 + b4*A4 + b2*A2 + b0*I
    W = np.dot(A6, W1) + W2
    U = np.dot(A, W)
    V = np.dot(A6, Z1) + Z2
    Lw1 = b13*M6 + b11*M4 + b9*M2
    Lw2 = b7*M6 + b5*M4 + b3*M2
    Lz1 = b12*M6 + b10*M4 + b8*M2
    Lz2 = b6*M6 + b4*M4 + b2*M2
    Lw = np.dot(A6, Lz1) + np.dot(M6, W1) + Lw2
    Lu = np.dot(A, Lw) + np.dot(E, W)
    Lv = np.dot(A6, Lz1) + np.dot(M6, Z1) + Lz2
    R, L = alg64_26_27(U, V, Lu, Lv)
    for k in range(1, s+1):
        L = np.dot(R, L) + np.dot(L, R)
        R = np.linalg.matrix_power(R, 2)
    return R, L

def main():
    for f in (linalg.expm, alg61):
        for M, E in zip(g_test_matrices, g_test_matrices_expm):
            print M
            print E
            print f(M)
            print '---------------------'
        print '==================================='

if __name__ == '__main__':
    main()
