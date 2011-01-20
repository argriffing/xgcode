"""
Construct matrices related to a distance matrix.

The matrices constructed here
are used in the proof of a theorem about partially supplied graphs.

The following notation will be used.
vertex counts and leaf multiplicity:
    q: number of leaves
    p: number of internal vertices
    N: the multiplicity of leaf vertices
vertex ordering abbreviation:
    LF: vertex ordering with leaf vertices first
    IF: vertex ordering with internal vertices first
matrix contents abbreviation:
    DO: original distance matrix (symmetric order p + q)
    DM: the expanded matrix (symmetric order p + Nq)
    DN: block of centered finitely expanded distance (symmetric order p + q)
    DI: block of centered infinitely expanded distance (symmetric order p + q)
matrix sub-block abbreviation:
    Q: leaf vertex principal submatrix (symmetric order q)
    P: internal vertex principal submatrix (symmetric order p)
    X: off-diagonal upper right submatrix (pxq or qxp)
"""

import unittest
import numpy as np

import NewickIO
import MatrixUtil
import EigUtil
import FelTree
import iterutils
import const

g_tree_string = const.read('20100730g').rstrip()


############################################################
# block matrix helper functions

def hrep(M, N):
    """
    Horizontally stack N copies of M.
    @param M: a matrix
    @param N: stack this many copies
    """
    return np.hstack([M]*N)

def vrep(M, N):
    """
    Vertically stack N copies of M.
    @param M: a matrix
    @param N: stack this many copies
    """
    return np.vstack([M]*N)

def get_corners(M, a, b):
    """
    Returns four matrices.
    The first is the top left axa matrix.
    The second is the top right axb matrix.
    The third is the bottom left bxa matrix.
    The fourth is the bottom right bxb matrix.
    @param M: a matrix
    @param a: the order of the first principal submatrix
    @param b: the order of the last principal submatrix
    @return: four matrices
    """
    return M[:a, :a], M[:a, -b:], M[-b:, :a], M[-b:, -b:]

def assemble_corners(A, B, C, D):
    """
    @param A: top left corner
    @param B: top right corner
    @param C: bottom left corner
    @param D: bottom right corner
    @return: an assembled matrix
    """
    try:
        return np.vstack([np.hstack([A, B]), np.hstack([C, D])])
    except ValueError as e:
        arr = [A.shape, B.shape, C.shape, D.shape]
        msg = ', '.join(str(x) for x in arr)
        raise ValueError(msg)

def tree_to_leaf_first_ids(tree):
    ordered_ids = []
    ordered_ids.extend(id(node) for node in tree.gen_tips())
    ordered_ids.extend(id(node) for node in tree.gen_internal_nodes())
    return ordered_ids


############################################################
# structured matrix definitions

class LFDO:
    def __init__(self, M, p, q):
        self.M = M
        self.p = p
        self.q = q

class LFDN:
    def __init__(self, M, p, q, N):
        self.M = M
        self.p = p
        self.q = q
        self.N = N

class LFDI:
    def __init__(self, M, p, q):
        self.M = M
        self.p = p
        self.q = q

class LFDM:
    def __init__(self, M, p, q, N):
        self.M = M
        self.p = p
        self.q = q
        self.N = N


############################################################
# structured matrix initialisation and conversion functions

def tree_string_to_LFDO(tree_string):
    tree = NewickIO.parse(tree_string, FelTree.NewickTree)
    nvertices = len(list(tree.preorder()))
    nleaves = len(list(tree.gen_tips()))
    ordered_ids = tree_to_leaf_first_ids(tree)
    D = np.array(tree.get_partial_distance_matrix(ordered_ids))
    q = nleaves
    p = nvertices - nleaves
    return LFDO(D, p, q)

def LFDO_to_LFDM(lfdo, N):
    D, p, q = lfdo.M, lfdo.p, lfdo.q
    Q, X, XT, P = get_corners(D, q, p)
    QX = vrep(hrep(Q, N), N)
    XX = vrep(X, N)
    PX = P
    M = assemble_corners(QX, XX, XX.T, PX)
    return LFDM(M, p, q, N)

def LFDM_to_LFDN(lfdm):
    M, p, q, N = lfdm.M, lfdm.p, lfdm.q, lfdm.N
    HMH = MatrixUtil.double_centered(M)
    v = p + q
    return LFDN(HMH[-v:, -v:], p, q, N)

def LFDO_to_LFDN(lfdo, N):
    D, p, q = lfdo.M, lfdo.p, lfdo.q
    weighted_colsums = N*D[:q].sum(axis=0) + D[q:].sum(axis=0)
    Q, X, XT, P = get_corners(D, q, p)
    weighted_grand_sum = N*N*np.sum(Q) + 2*N*np.sum(X) + np.sum(P)
    weighted_colmeans = weighted_colsums / float(p+N*q)
    weighted_grand_mean = weighted_grand_sum / float(N*N*q*q + 2*N*p*q + p*p)
    M = D.copy()
    M -= weighted_colmeans
    M = M.T
    M -= weighted_colmeans
    M += weighted_grand_mean
    return LFDN(M, p, q, N)

def LFDO_to_LFDI(lfdo):
    D, p, q = lfdo.M, lfdo.p, lfdo.q
    weighted_colmeans = D[:q].mean(axis=0)
    Q, X, XT, P = get_corners(D, q, p)
    weighted_grand_mean = np.mean(Q)
    M = D.copy()
    M -= weighted_colmeans
    M = M.T
    M -= weighted_colmeans
    M += weighted_grand_mean
    return LFDI(M, p, q)


############################################################
# test invariants
# test shortcuts

class ProofDecorationTest(unittest.TestCase):

    def test_lfdn_shortcut(self):
        """
        Check two ways of constructing blocks of a matrix.
        The matrix in question is the centered finitely extended matrix.
        One way is to construct the whole matrix and then take blocks.
        The other way is to construct the blocks more directly.
        """
        lfdo = tree_string_to_LFDO(g_tree_string)
        for N in (1, 5, 10):
            lfdn_a = LFDO_to_LFDN(lfdo, N)
            lfdm = LFDO_to_LFDM(lfdo, N)
            lfdn_b = LFDM_to_LFDN(lfdm)
            self.assertTrue(np.allclose(lfdn_a.M, lfdn_b.M))

    def test_degenerate_lfdn(self):
        """
        Make sure that replication with a factor of 1 does not do anything.
        """
        lfdo = tree_string_to_LFDO(g_tree_string)
        N = 1
        lfdn = LFDO_to_LFDN(lfdo, N)
        HDH = MatrixUtil.double_centered(lfdo.M)
        self.assertTrue(np.allclose(lfdn.M, HDH))

    def test_lfdi_approximation(self):
        """
        As N increases, the approximation should become closer.
        More precisely, as N becomes large,
        multiplying N by ten should
        add one decimal place of accuracy to the approximation.
        Where the accuracy of the approximation is taken
        to be the frobenius norm of the error matrix.
        """
        lfdo = tree_string_to_LFDO(g_tree_string)
        lfdi = LFDO_to_LFDI(lfdo)
        # For these values of N,
        # the error for N should be more than 9 times the error for 10N.
        # When N is very large,
        # the error for N should approach 10 times the error for 10N.
        Ns = (10, 100, 1000, 10000)
        lfdns = [LFDO_to_LFDN(lfdo, N) for N in Ns]
        error_norms = [np.linalg.norm(lfdi.M - lfdn.M) for lfdn in lfdns]
        for ea, eb in iterutils.pairwise(error_norms):
            # ea should be more than nine times as bad as eb
            self.assertTrue(ea / eb > 9)

    def test_lfdi_submatrix(self):
        """
        Test a principal submatrix of the LFDI matrix.
        The principal submatrix of the LFDI which corresponds to leaf vertices
        should equal the centered leaf vertex distance matrix.
        """
        lfdo = tree_string_to_LFDO(g_tree_string)
        q = lfdo.q
        DQ = lfdo.M[:q, :q]
        HDQH = MatrixUtil.double_centered(DQ)
        lfdi = LFDO_to_LFDI(lfdo)
        self.assertTrue(np.allclose(lfdi.M[:q, :q], HDQH))


if __name__ == '__main__':
    unittest.main()
