"""
These are matrix operations that depend primarily on the numpy package.
Matrices in this module are 2d numpy arrays, not numpy Matrix objects.
Vectors are 1d numpy arrays, not 1xN or Nx1 2d structures.
"""

import numpy as np

import MatrixUtil
from MatrixUtil import MatrixError

def bott_duffin(M, v):
    """
    Compute a constrained generalized inverse.
    Specifically, this is the Bott-Duffin inverse of M
    constrained to the orthogonal complement of v.
    This function assumes that v has rank 1,
    although Bott-Duffin inverses are also defined
    for inverses constrained to orthogonal complements
    of higher dimensional subspaces.
    Maybe this could be a separate python function
    where v is replaced by a shape-2 numpy array.
    @param M: a matrix
    @param v: a vector
    @return: the constrained generalized inverse of M
    """
    # check the shapes of the input matrix and vector
    MatrixUtil.assert_1d(v)
    n = len(v)
    if M.shape != (n, n):
        raise ValueError('M and v have incompatible shapes')
    # check that v is nonzero
    v_dot_v = np.inner(v, v)
    if not v_dot_v:
        raise ValueError('expected nonzero v')
    # compute the orthogonal projection onto v
    P = np.outer(v, v) / v_dot_v
    # compute the orthogonal projection onto the orthogonal complement of v
    I = np.eye(n)
    C = I - P
    # compute the constrained generalized inverse
    B = np.dot(C, np.linalg.inv(np.dot(M, C) + P))
    return B

def bott_duffin_const(M):
    """
    Compute a constrained generalized inverse.
    Specifically, this is the Bott-Duffin inverse of M
    constrained to the orthogonal complement of the constant vector.
    """
    MatrixUtil.assert_square(M)
    n = len(M)
    e = np.ones(n)
    return bott_duffin(M, e)
