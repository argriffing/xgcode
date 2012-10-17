#
# Authors: Travis Oliphant, March 2002
#          Anthony Scopatz, August 2012 (Sparse Updates)
#          Jake Vanderplas, August 2012 (Sparse Updates)
#
# Also includes expm test data from
# http://people.sc.fsu.edu/~jburkardt/m_src/test_matrix_exponential
#

import math

import numpy as np
from scipy import linalg

# Use the correct eps when understanding Matlab notation.
# http://www.scipy.org/NumPy_for_Matlab_Users
matlab_eps = np.spacing(1)

# More constants for the expm test matrices.
exp16 = math.exp(16)
exp4 = math.exp(4)


class ExpmCase:
    def __init__(self, burkardt_index, description, M, expm_M):
        self.burkardt_index = burkardt_index
        self.description = description
        self.M = M
        self.expm_M = expm_M

g_burkardt_expm_test_data = [

        ExpmCase(
            1,
            'This matrix is diagonal.\n'
            'The calculation of the matrix exponential is simple.',
            np.array([
                [1, 0],
                [0, 2],
                ], dtype=float),
            np.array([
                [2.718281828459046, 0],
                [0, 7.389056098930650],
                ], dtype=float),
            ),

        ExpmCase(
            2,
            'This matrix is symmetric.\n'
            'The calculation of the matrix exponential is straightforward.',
            np.array([
                [1, 3],
                [3, 2],
                ], dtype=float),
            np.array([
                [39.322809708033859, 46.166301438885753],
                [46.166301438885768, 54.711576854329110],
                ], dtype=float),
            ),

        ExpmCase(
            3,
            'This example is due to Laub.\n'
            'This matrix is ill-suited for the Taylor series approach.\n'
            'As powers of A are computed, the entries blow up too quickly.',
            np.array([
                [0, 1],
                [-39, -40],
                ], dtype=float),
            np.array([
                [0, 2.718281828459046],
                [1.154822e-17, 2.718281828459046],
                ], dtype=float),
            ),

        ExpmCase(
            4,
            'This example is due to Moler and Van Loan.\n'
            'The example will cause problems '
            'for the series summation approach,\n'
            'as well as for diagonal Pade approximations.',
            np.array([
                [-49, 24],
                [-64, 31],
                ], dtype=float),
            np.array([
                [-0.735759, 0.551819],
                [-1.471518, 1.103638],
                ], dtype=float),
            ),

        ExpmCase(
            5,
            'This example is due to Moler and Van Loan.\n'
            'This matrix is strictly upper triangular\n'
            'All powers of A are zero beyond some (low) limit.\n'
            'This example will cause problems for Pade approximations.',
            np.array([
                [0, 6, 0, 0],
                [0, 0, 6, 0],
                [0, 0, 0, 6],
                [0, 0, 0, 0],
                ], dtype=float),
            np.array([
                [1, 6, 18, 36],
                [0, 1, 6, 18],
                [0, 0, 1, 6],
                [0, 0, 0, 1],
                ], dtype=float),
            ),

        ExpmCase(
            6,
            'This example is due to Moler and Van Loan.\n'
            'This matrix does not have a complete set of eigenvectors.\n'
            'That means the eigenvector approach will fail.',
            np.array([
                [1, 1],
                [0, 1],
                ], dtype=float),
            np.array([
                [2.718281828459046, 2.718281828459046],
                [0, 2.718281828459046],
                ], dtype=float),
            ),

        ExpmCase(
            7,
            'This example is due to Moler and Van Loan.\n'
            'This matrix is very close to example 5.\n'
            'Mathematically, it has a complete set of eigenvectors.\n'
            'Numerically, however, the calculation will be suspect.',
            np.array([
                [1 + matlab_eps, 1],
                [0, 1 - matlab_eps],
                ], dtype=float),
            np.array([
                [2.718309, 2.718282],
                [0, 2.718255],
                ], dtype=float),
            ),

        ExpmCase(
            8,
            'This matrix was an example in Wikipedia.',
            np.array([
                [21, 17, 6],
                [-5, -1, -6],
                [4, 4, 16],
                ], dtype=float),
            np.array([
                [13*exp16 - exp4, 13*exp16 - 5*exp4,  2*exp16 - 2*exp4],
                [-9*exp16 + exp4, -9*exp16 + 5*exp4, -2*exp16 + 2*exp4],
                [16*exp16,        16*exp16,           4*exp16         ],
                ], dtype=float) * 0.25,
            ),

        ExpmCase(
            9,
            'This matrix is due to the NAG Library.\n'
            'It is an example for function F01ECF.',
            np.array([
                [1, 2, 2, 2],
                [3, 1, 1, 2],
                [3, 2, 1, 2],
                [3, 3, 3, 1],
                ], dtype=float),
            np.array([
                [740.7038, 610.8500, 542.2743, 549.1753],
                [731.2510, 603.5524, 535.0884, 542.2743],
                [823.7630, 679.4257, 603.5524, 610.8500],
                [998.4355, 823.7630, 731.2510, 740.7038],
                ], dtype=float),
            ),

        ExpmCase(
            10,
            'This is Ward\'s example #1.\n'
            'It is defective and nonderogatory.\n'
            'The eigenvalues are 3, 3 and 6.',
            np.array([
                [4, 2, 0],
                [1, 4, 1],
                [1, 1, 4],
                ], dtype=float),
            np.array([
                [147.8666224463699, 183.7651386463682, 71.79703239999647],
                [127.7810855231823, 183.7651386463682, 91.88256932318415],
                [127.7810855231824, 163.6796017231806, 111.9681062463718],
                ], dtype=float),
            ),

        ExpmCase(
            11,
            'This is Ward\'s example #2.\n'
            'It is a symmetric matrix.\n'
            'The eigenvalues are 20, 30, 40.',
            np.array([
                [29.87942128909879, 0.7815750847907159, -2.289519314033932],
                [0.7815750847907159, 25.72656945571064, 8.680737820540137],
                [-2.289519314033932, 8.680737820540137, 34.39400925519054],
                ], dtype=float),
            np.array([
                 [
                     5.496313853692378E+15,
                     -1.823188097200898E+16,
                     -3.047577080858001E+16],
                 [
                    -1.823188097200899E+16,
                    6.060522870222108E+16,
                    1.012918429302482E+17],
                 [
                    -3.047577080858001E+16,
                    1.012918429302482E+17,
                    1.692944112408493E+17],
                ], dtype=float),
            ),

        ExpmCase(
            12,
            'This is Ward\'s example #3.\n'
            'Ward\'s algorithm has difficulty estimating the accuracy\n'
            'of its results.  The eigenvalues are -1, -2, -20.',
            np.array([
                [-131, 19, 18],
                [-390, 56, 54],
                [-387, 57, 52],
                ], dtype=float),
            np.array([
                [-1.509644158793135, 0.3678794391096522, 0.1353352811751005],
                [-5.632570799891469, 1.471517758499875, 0.4060058435250609],
                [-4.934938326088363, 1.103638317328798, 0.5413411267617766],
                ], dtype=float),
            ),
        ]


# Skip this one because it is parameterized.
"""
        BurkardtExpmTest(
            13,
            'This is Ward\'s example #4.\n'
            'This is a version of the Forsythe matrix.\n'
            'The eigenvector problem is badly conditioned.\n'
            'Ward\'s algorithm has difficulty estimating the accuracy\n'
            'of its results for this problem.'
            foo,
            bar,
            ),
        ]
"""


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
    Reference is
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

