"""
Try to implement Algorithm 6.4 of a paper by Awad Al-Mohy and Nick Higham.

Computing the Frechet derivative of the matrix exponential,
with an application to condition number estimation.
The magic constants in this module are tuned
for double precision floating point computations.
"""

import numpy as np
from scipy import linalg

# Coefficients of degree 13 Pade approximant.
# From Algorithm 2.3 of Higham 2005,
# The Scaling and Squaring Method for the Matrix Exponential Revisited.
# b0 through b13
btable = [
        64752532480000, 32382376266240000, 7771770303897600,
        1187353796428800, 129060195264000, 10559470521600,
        670442572800, 33522128640, 1323241920,
        40840800, 960960, 16380, 182, 1,
        ]
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

def alg611(X, m):
    """
    This is a naive transcription of definition (6.11).
    It is a decomposition of a sum into odd and even parts.
    """
    U = np.zeros_like(X)
    V = np.zeros_like(X)
    for k in range(0, (m-1)//2 + 1):
        X2k = linalg.matrix_power(X, 2*k)
        U += btable[2*k + 1] * X2k
        V += btable[2*k] * X2k
    return np.dot(X, U), V

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
        return np.dot(alg34(k-1), A) + np.dot(linalg.matrix_power(A, k-1), E)

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
        Lub += btable[2*k + 1] * linalg.matrix_power(A, 2*k)
    Lu = np.dot(A, Lua) + np.dot(E, Lub)
    return Lu, Lv

def alg64(A, E):
    """
    Scaling and squaring algorithm for exponential and Frechet derivative.
    Given A, E, this algorithm computes R and L by scaling and squaring.
    The matrices A and E are square nxn with complex entries.
    Matrix R is expm(A), and matrix L is Lexp(A, E) in their notation.
    """
    for m in (3, 5, 7, 9):
        if linalg.norm(A, 1) < table61_l[m]:
            # Evaluate U = um(A) and V = vm(A), using (6.11).
            U, V = alg611(A, m)
            # Evaluate Lu = Lum(A, E) and Lv = Lvm(A, E),
            # using (6.12) and (6.13).
            Lu, Lv = alg612(A, E, m)
            #XXX


def main():
    pass

if __name__ == '__main__':
    main()
