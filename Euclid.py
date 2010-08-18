"""
This module is about spherical Euclidean distance matrices.
A Euclidean distance matrix is defined here as in the
paper "On Euclidean distance matrices" by Balaji and Bapat.
That is, it is a symmetric matrix of squared Euclidean distances.
This paper's definition of a Laplacian is also used here.
Where matrices are used in this module, they are implemented as numpy arrays.
"""

import unittest
import random

import numpy as np

import MatrixUtil
import SchurAlgebra

# presumably this is from some tree
g_D_a = np.array([
    [0, 3, 3, 3, 3],
    [3, 0, 2, 4, 4],
    [3, 2, 0, 4, 4],
    [3, 4, 4, 0, 2],
    [3, 4, 4, 2, 0]], dtype=float)

# this is from the ((0:1, 1:2):3, 2:4, 3:5); tree
g_D_b = np.array([
    [0,  3, 8,  9],
    [3,  0, 9, 10],
    [8,  9, 0,  9],
    [9, 10, 9,  0]], dtype=float)

# this is from the ((0:1, 1:2):3, 2:4, 3:5); tree and includes internal nodes
g_D_c = np.array([
    [0,  3, 8,  9, 1, 4],
    [3,  0, 9, 10, 2, 5],
    [8,  9, 0,  9, 7, 4],
    [9, 10, 9,  0, 8, 5],
    [1,  2, 7,  8, 0, 3],
    [4,  5, 4,  5, 3, 0]], dtype=float)

def edm_to_adjacency(D):
    """
    @param D: a euclidean distance matrix
    @return: an adjacency matrix
    """
    return laplacian_to_adjacency(edm_to_laplacian(D))

def laplacian_to_adjacency(L):
    """
    @param L: a laplacian matrix
    @return: an adjacency matrix
    """
    MatrixUtil.assert_square(L)
    return np.diag(np.diag(L)) - L

def adjacency_to_laplacian(A):
    """
    @param A: an adjacency matrix
    @return: a laplacian matrix
    """
    MatrixUtil.assert_square(A)
    return np.diag(np.sum(A, 0)) - A

def cov_to_edm(S):
    """
    @param S: a covariance matrix
    @return: a Euclidean distance matrix
    """
    MatrixUtil.assert_square(S)
    n = len(S)
    d = np.diag(S)
    e = np.ones_like(d)
    D = np.outer(d, e) + np.outer(e, d) - 2*S
    return D

def dccov_to_edm(HSH):
    """
    @param HSH: a double centered covariance matrix
    @return: a Euclidean distance matrix
    """
    MatrixUtil.assert_square(HSH)
    return cov_to_edm(HSH)

def laplacian_to_dccov(L):
    """
    @param L: a laplacian matrix
    @return: a double centered covariance matrix
    """
    MatrixUtil.assert_square(L)
    M = np.ones_like(L) / float(len(L))
    # This should be the same but perhaps not as numerically stable:
    # HSH = np.linalg.pinv(L)
    HSH = np.linalg.pinv(L - M) + M
    return HSH

def dccov_to_laplacian(HSH):
    """
    This function stably pseudoinverts a double centered matrix.
    @param HSH: a double centered covariance matrix
    @return: a laplacian matrix
    """
    MatrixUtil.assert_square(HSH)
    M = np.ones_like(HSH) / float(len(HSH))
    # This should be the same but perhaps not as numerically stable:
    # L = np.linalg.pinv(HSH)
    L = np.linalg.pinv(HSH - M) + M
    return L

def edm_to_dccov(D):
    """
    @param D: a Euclidean distance matrix
    @return: a double centered covariance matrix
    """
    MatrixUtil.assert_square(D)
    return -(0.5)*MatrixUtil.double_centered(D)

def laplacian_to_edm(L):
    """
    @param L: a laplacian matrix
    @return: a Euclidean distance matrix
    """
    MatrixUtil.assert_square(L)
    return dccov_to_edm(laplacian_to_dccov(L))

def edm_to_laplacian(D):
    """
    @param D: a Euclidean distance matrix
    @return: a laplacian matrix
    """
    MatrixUtil.assert_square(D)
    D_pinv = np.linalg.pinv(D)
    a = np.sum(D_pinv, 0)
    denom = np.sum(D_pinv)
    if not denom:
        raise ValueError('the EDM is not spherical')
    L = -2.0*(D_pinv - np.outer(a,a)/denom)
    return L

def edm_to_q(D):
    """
    @param D: a treelike distance matrix
    @return: the neighbor joining Q matrix
    """
    MatrixUtil.assert_square(D)
    n = len(D)
    r = np.sum(D, 0)
    e = np.ones_like(r)
    Q = (n-2)*D - np.outer(e, r) - np.outer(r, e)
    return Q

def dccov_to_q(HSH):
    """
    Is there a shortcut here?
    @param HSH: a doubly centered covariance matrix
    @return: a neighbor joining Q matrix
    """
    MatrixUtil.assert_square(HSH)
    n = len(HSH)
    d = np.diag(HSH)
    e = np.ones_like(d)
    D = np.outer(d, e) + np.outer(e, d) - 2*HSH
    r = np.sum(D, 0)
    e = np.ones_like(r)
    Q = (n-2)*D - np.outer(e, r) - np.outer(r, e)
    return Q

def q_to_dccov(Q):
    """
    @param Q: a neighbor joining Q matrix
    @return: something like a doubly centered covariance matrix
    """
    return MatrixUtil.double_centered(q_to_cov(Q))

def q_to_cov(Q):
    """
    @param Q: a neighbor joining Q matrix
    @return: something like a covariance matrix
    """
    MatrixUtil.assert_square(Q)
    n = len(Q)
    S = -Q/(2*(n-2))
    return S

def q_to_edm(Q):
    """
    @param Q: a neighbor joining Q matrix
    @return: A Euclidean distance matrix
    """
    return cov_to_edm(q_to_cov(Q))

def dccov_to_points(HSH, naxes=None):
    """
    Get points where the first coordinate is on the principal axis.
    @param HSH: a centered covariance matrix or a Gower matrix
    @param naxes: use this many axes, defaulting to N-1
    @return: a matrix whose rows define N points in (naxes)-space
    """
    if naxes is None:
        naxes = len(HSH) - 1
    W, VT = np.linalg.eigh(HSH)
    V = VT.T.tolist()
    vectors = [np.array(v)*w for w, v in list(reversed(sorted(zip(np.sqrt(W), V))))[:naxes]]
    X = np.array(zip(*vectors))
    return X

def edm_to_points(D, naxes=None):
    """
    Get points where the first coordinate is on the principal axis.
    @param D: a matrix of squared Euclidean distances
    @param naxes: use this many axes, defaulting to N-1
    @return: a matrix whose rows define N points in (naxes)-space
    """
    HSH = edm_to_dccov(D)
    return dccov_to_points(HSH, naxes)

def edm_to_weighted_cross_product(D, m):
    """
    Get the cross product matrix.
    This uses the terminology of Herve Abdi 2007.
    @param D: a matrix of squared Euclidean distances
    @param m: a vector defining the mass distribution
    @return: a cross product matrix
    """
    n = len(m)
    if D.shape != (n, n):
        raise ValueError('D should be a square matrix conformant to m')
    if any(x < 0 for x in m):
        raise ValueError('each element in m should be nonnegative')
    if not np.allclose(np.sum(m), 1):
        raise ValueError('the masses should sum to one')
    # construct the centering matrix
    I = np.eye(n, dtype=float)
    E = I - np.outer(np.ones(n, dtype=float), m)
    # get the cross product matrix
    S = (-0.5)*np.dot(E, np.dot(D, E.T))
    return S

def edm_to_weighted_points(D, m):
    """
    @param D: a distance matrix
    @param m: a mass vector
    @return: a matrix whose rows are embedded points
    """
    # Here is an excerpt of the email I sent to Eric on 20091007
    # 
    # Let D be a matrix of squared pairwise Euclidean distances.  Let m be a
    # column vector of nonnegative weights such that the weights sum to 1.
    # Now we create the weighted centering matrix E = I - em'.  Note that if
    # the weights are uniform then this is the standard centering matrix H.
    # Next define Y such that the positive semidefinite matrix YY' =
    # -(1/2)EDE'.  Note that if D is NxN there is a Nx(N-1) Y matrix that
    # satisfies this equality.  Now find the singular value decomposition USV'
    # = MY, and compute Z = YV.  Rows of Z are the points in N-1 dimensional
    # space where the basis vectors are the (weighted) MDS axes.
    # 
    # define the total number of vertices
    nvertices = len(m)
    # get the cross product (S in the Abdi notation)
    cross_product = edm_to_weighted_cross_product(D, m)
    # get rows of embedded points whose origin is defined by the mass vector
    W, VT = np.linalg.eigh(cross_product)
    V = VT.T.tolist()
    vectors = [np.array(v)*w for w, v in list(reversed(sorted(zip(np.sqrt(W), V))))[:nvertices-1]]
    XT = np.array(vectors)
    Y = XT.T
    # do this projection
    M = np.diag(np.sqrt(m))
    U, S, VT = np.linalg.svd(np.dot(M, Y))
    Z = np.dot(Y, VT.T)
    return Z


class TestEuclid(unittest.TestCase):

    def nonzero_reciprocal(self, M):
        n = len(M)
        R = M.copy()
        for i in range(n):
            for j in range(n):
                v = R[i][j]
                R[i][j] = (0.0 if not v else 1.0 / v)
        return R

    def test_edm_to_laplacian_and_back(self):
        D = np.array([
            [0.0, 3.0, 3.0, 3.0, 3.0],
            [3.0, 0.0, 2.0, 4.0, 4.0],
            [3.0, 2.0, 0.0, 4.0, 4.0],
            [3.0, 4.0, 4.0, 0.0, 2.0],
            [3.0, 4.0, 4.0, 2.0, 0.0]])
        L = edm_to_laplacian(D)
        D_estimated = laplacian_to_edm(L)
        self.assertTrue(np.allclose(D, D_estimated))

    def test_adjacency_to_laplacian(self):
        A = np.array([
            [0, 0, 0, 1],
            [0, 0, 0, 2],
            [0, 0, 0, 3],
            [1, 2, 3, 0]])
        L = np.array([
            [1, 0, 0, -1],
            [0, 2, 0, -2],
            [0, 0, 3, -3],
            [-1, -2, -3, 6]])
        expected = L
        observed = adjacency_to_laplacian(A)
        self.assertTrue(np.allclose(observed, expected))

    def test_laplacian_to_adjacency(self):
        A = np.array([
            [0, 0, 0, 1],
            [0, 0, 0, 2],
            [0, 0, 0, 3],
            [1, 2, 3, 0]])
        L = np.array([
            [1, 0, 0, -1],
            [0, 2, 0, -2],
            [0, 0, 3, -3],
            [-1, -2, -3, 6]])
        expected = A
        observed = laplacian_to_adjacency(L)
        self.assertTrue(np.allclose(observed, expected))

    def test_commutativity(self):
        """
        Schur complementation and merging can be done in either order.
        """
        reciprocal_adjacency_big = np.array([
                [0, 0, 0, 0, 0, 0, 0, 2, 0, 0],
                [0, 0, 0, 0, 0, 0, 2, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 9, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 3, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 2],
                [0, 2, 9, 0, 0, 0, 0, 4, 0, 0],
                [2, 0, 0, 0, 0, 0, 4, 0, 0, 1],
                [0, 0, 0, 1, 3, 0, 0, 0, 0, 7],
                [0, 0, 0, 0, 0, 2, 0, 1, 7, 0]], dtype=float)
        A_big = self.nonzero_reciprocal(reciprocal_adjacency_big)
        L_big = adjacency_to_laplacian(A_big)
        # define the pruned branch length
        p = 101.0 / 39.0
        reciprocal_adjacency_small = np.array([
                [0, 0, 0, 0, 0, 2],
                [0, 0, 0, 0, 2, 0],
                [0, 0, 0, 0, 9, 0],
                [0, 0, 0, 0, 0, p],
                [0, 2, 9, 0, 0, 4],
                [2, 0, 0, p, 4, 0]])
        A_small = self.nonzero_reciprocal(reciprocal_adjacency_small)
        L_small = adjacency_to_laplacian(A_small)
        # get the small matrix in terms of the big matrix by schur complementation followed by merging
        reconstructed_small_a = SchurAlgebra.mmerge(SchurAlgebra.mschur(L_big, set([8, 9])), set([3, 4, 5]))
        self.assertTrue(np.allclose(L_small, reconstructed_small_a))
        # get the small matrix in terms of the big matrix by merging followed by schur complementation
        reconstructed_small_b = SchurAlgebra.mschur(SchurAlgebra.mmerge(L_big, set([3, 4, 5])), set([6, 7]))
        self.assertTrue(np.allclose(L_small, reconstructed_small_b))
        # get the laplacian associated with a 4x4 distance matrix in multiple ways
        first_result = SchurAlgebra.mmerge(SchurAlgebra.mschur(L_big, set([6, 7, 8, 9])), set([3, 4, 5]))
        second_result = SchurAlgebra.mschur(L_small, set([4, 5]))
        self.assertTrue(np.allclose(first_result, second_result))
    
    def test_edm_to_q_and_back_long(self):
        """
        Go from D to Q to HSH to D.
        """
        D = np.array([
            [0.0, 3.0, 3.0, 3.0, 3.0],
            [3.0, 0.0, 2.0, 4.0, 4.0],
            [3.0, 2.0, 0.0, 4.0, 4.0],
            [3.0, 4.0, 4.0, 0.0, 2.0],
            [3.0, 4.0, 4.0, 2.0, 0.0]])
        Q = edm_to_q(D)
        HSH = q_to_dccov(Q)
        D_prime = dccov_to_edm(HSH)
        self.assertTrue(np.allclose(D, D_prime))
    
    def test_edm_to_q_and_back_short(self):
        """
        Go from D to Q to D.
        """
        D = np.array([
            [0.0, 3.0, 3.0, 3.0, 3.0],
            [3.0, 0.0, 2.0, 4.0, 4.0],
            [3.0, 2.0, 0.0, 4.0, 4.0],
            [3.0, 4.0, 4.0, 0.0, 2.0],
            [3.0, 4.0, 4.0, 2.0, 0.0]])
        Q = edm_to_q(D)
        D_prime = q_to_edm(Q)
        self.assertTrue(np.allclose(D, D_prime))

    def test_edm_to_points(self):
        """
        Test the construction of points from a matrix of squared distances.
        Assert that they give back the distances from which they were constructed.
        """
        # check a couple of distance matrices
        for D in (g_D_a, g_D_b):
            X = edm_to_points(D)
            # construct the squared distances directly from the points
            D_prime = np.array([[np.dot(pb-pa, pb-pa) for pa in X] for pb in X])
            self.assertTrue(np.allclose(D, D_prime))
            # check another property of the matrix of embedded points
            D_prime = dccov_to_edm(np.dot(X, X.T))
            self.assertTrue(np.allclose(D, D_prime))

    def test_edm_to_weighted_points(self):
        """
        Test the construction of points from a matrix of pairwise squared distances among weighted points.
        """
        # use a distance matrix from an asymmetric tree and include internal nodes
        D = g_D_c
        n = len(D)
        # try some different hardcoded mass distributions
        self.assertTrue(n == 6)
        m_uniform = np.ones(n) / float(n)
        m_weighted_positive = np.array([.1, .1, .1, .1, .3, .3])
        m_weighted_nonnegative = np.array([.25, .25, .25, .25, 0.0, 0.0])
        # assert some invariants
        for m in (m_uniform, m_weighted_positive, m_weighted_nonnegative):
            # get the points
            X = edm_to_weighted_points(D, m)
            # check one way of getting the distance matrix back from the points
            D_prime = dccov_to_edm(np.dot(X, X.T))
            msg = (D, D_prime)
            self.assertTrue(np.allclose(D, D_prime), msg)
            # Check one way of getting the distance matrix back from the cross product matrix.
            # This is equation (4) of the 2007 Abdi chapter.
            S = edm_to_weighted_cross_product(D, m)
            D_prime = dccov_to_edm(S)
            msg = (D, D_prime)
            self.assertTrue(np.allclose(D, D_prime), msg)
            # compute the pairwise distances of these points
            D_prime = np.array([[np.dot(pb-pa, pb-pa) for pa in X] for pb in X])
            msg = (D, D_prime)
            self.assertTrue(np.allclose(D, D_prime), msg)


if __name__ == '__main__':
    np.set_printoptions(linewidth=200)
    unittest.main()
