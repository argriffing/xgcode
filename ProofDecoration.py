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

from SnippetUtil import HandlingError
import SnippetUtil
import Form
import FormOut
import NewickIO
import Euclid
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
        lfdo = tree_string_to_LFDO(g_tree_string)
        for N in (1, 5, 10):
            lfdn_a = LFDO_to_LFDN(lfdo, N)
            lfdm = LFDO_to_LFDM(lfdo, N)
            lfdn_b = LFDM_to_LFDN(lfdm)
            self.assertTrue(np.allclose(lfdn_a.M, lfdn_b.M))

    def test_degenerate_lfdn(self):
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





############################################################
# FIXME old obsolete stuff that should be deleted


class DistanceModel:
    """
    This is a pure base class.
    """

    def __init__(self, tree_string):
        tree = NewickIO.parse(tree_string, FelTree.NewickTree)
        nvertices = len(list(tree.preorder()))
        nleaves = len(list(tree.gen_tips()))
        ordered_ids = self.tree_to_ordered_ids(tree)
        self.D = np.array(tree.get_partial_distance_matrix(ordered_ids))
        self.q = nleaves
        self.p = nvertices - nleaves

    def get_Dpq(self):
        return self.D, self.p, self.q

def LF_get_big_D(D, q, N):
    """
    LF means leaf first.
    Note that D should be returned for N=1.
    @param D: a distance matrix with leaves first
    @param q: the number of leaves
    @param N: expansion factor
    """
    if N < 2:
        return self.D
    p = len(D) - q
    DQ, DX, DXT, DP = get_corners(D, q, p)
    DQ_exp = vrep(hrep(DQ, N), N)
    DX_exp = vrep(DX, N)
    DP_exp = DP
    return assemble_corners(DQ_exp, DX_exp, DX_exp.T, DP_exp)

class LeafFirstModel(DistanceModel):

    def tree_to_ordered_ids(self, tree):
        ordered_ids = []
        ordered_ids.extend(id(node) for node in tree.gen_tips())
        ordered_ids.extend(id(node) for node in tree.gen_internal_nodes())
        return ordered_ids

    def get_DQ(self):
        D, p, q = self.get_Dpq()
        return D[:q, :q]

    def get_DP(self):
        D, p, q = self.get_Dpq()
        return D[q:, q:]

    def get_DX(self):
        D, p, q = self.get_Dpq()
        return D[:q, q:]

    def get_expanded_matrix(self, N):
        """
        Note that D should be returned for N=1.
        @param N: the number of replicates
        @return: a distance matrix
        """
        if N < 2:
            return self.D
        DQ, DP, DX = self.get_DQ(), self.get_DP(), self.get_DX()
        DQ_exp = np.hstack([np.vstack([DQ]*N)]*N)
        DX_exp = np.vstack([DX]*N)
        DP_exp = DP
        return hstack([vstack([DQ_exp, DX_exp.T]),
            vstack([DX_exp, DP_exp])])

    def get_centered_expanded_matrix(N):
        DX = self.get_expanded_matrix(N)
        return MatrixUtil.double_centered(DX)

    def get_centered_expanded_DQ(N):
        """
        Get a block of the expanded and centered matrix.
        """
        p, q = self.p, self.q
        D = get_centered_expanded_matrix(N)
        return D[:q, :q]

    def get_centered_expanded_DX(N):
        """
        Get a block of the expanded and centered matrix.
        """
        p, q = self.p, self.q
        D = get_centered_expanded_matrix(N)
        return D[:q, :-p]


class InternalFirstModel(DistanceModel):

    def tree_to_ordered_ids(self, tree):
        ordered_ids = []
        ordered_ids.extend(id(node) for node in tree.gen_internal_nodes())
        ordered_ids.extend(id(node) for node in tree.gen_tips())
        return ordered_ids

    def get_DQ(self):
        D, p, q = self.get_Dpq()
        return D[p:, p:]

    def get_DP(self):
        D, p, q = self.get_Dpq()
        return D[:p, :p]

    def get_DX(self):
        D, p, q = self.get_Dpq()
        return D[:p, p:]


def get_grant_proposal_points(D, nleaves):
    """
    @return: rows are 2D points
    """
    B = D_to_B_third(D, nleaves)
    DQ = D_to_DQ_internal_first(D, nleaves)
    GQ = (-0.5) * MatrixUtil.double_centered(DQ)
    ws, vs = EigUtil.eigh(GQ)
    S = np.diag(ws)
    U = np.vstack(vs).T
    USUT = np.dot(np.dot(U, S), U.T)
    if not np.allclose(USUT, GQ):
        raise ValueError('eigenfail')
    S_sqrt = np.diag(np.sqrt(ws))
    S_sqrt_pinv = np.linalg.pinv(S_sqrt)
    # leaf-mds points
    X = np.dot(U, S_sqrt)
    # imputed internal points
    # the following line is as written in the proposal but it is wrong
    #W = np.dot(np.dot(S_sqrt_pinv, B), U)
    try:
        W = np.dot(np.dot(B, U), S_sqrt_pinv)
    except ValueError as e:
        arr = [
                B.shape,
                U.shape,
                S_sqrt_pinv.shape]
        msg = ', '.join(str(x) for x in arr)
        raise ValueError(msg)
    # put them together and get only the first coordinates
    full_points = np.vstack([W, X])
    points = full_points.T[:2].T
    return points

def get_index_edges(tree, ordered_ids):
    """
    Given a tree and some ordered ids, get edges defined on indices.
    @param tree: the tree object
    @param ordered_ids: the returned index pairs are for this sequence
    @return: a collection of index pairs defining edges
    """
    # map ids to indices
    id_to_index = dict((myid, index) for index, myid in enumerate(ordered_ids))
    # each edge in this set is a frozenset of two indices
    index_edges = set()
    for node in tree.preorder():
        index = id_to_index[id(node)]
        for neighbor in node.gen_neighbors():
            neighbor_index = id_to_index[id(neighbor)] 
            index_edges.add(frozenset([index, neighbor_index]))
    return index_edges

def get_animation_frame(
        image_format, physical_size, scale, index_edges, points):
    """
    This function is about drawing the tree.
    @param image_format: the image extension
    @param physical_size: the width and height of the image in pixels
    @param scale: a scaling factor
    @param index_edges: defines the connectivity of the tree
    @param points: an array of 2D points, the first few of which are leaves
    @return: the animation frame as an image as a string
    """
    # before we begin drawing we need to create the cairo surface and context
    cairo_helper = CairoUtil.CairoHelper(image_format)
    surface = cairo_helper.create_surface(physical_size[0], physical_size[1])
    context = cairo.Context(surface)
    # define some helper variables
    x0 = physical_size[0] / 2.0
    y0 = physical_size[1] / 2.0
    npoints = len(points)
    # draw an off-white background
    context.save()
    context.set_source_rgb(.9, .9, .9)
    context.paint()
    context.restore()
    # draw the axes which are always in the center of the image
    context.save()
    context.set_source_rgb(.9, .7, .7)
    context.move_to(x0, 0)
    context.line_to(x0, physical_size[1])
    context.stroke()
    context.move_to(0, y0)
    context.line_to(physical_size[0], y0)
    context.stroke()
    context.restore()
    # draw the edges
    context.save()
    context.set_source_rgb(.8, .8, .8)
    for edge in index_edges:
        ai, bi = tuple(edge)
        ax, ay = points[ai].tolist()
        bx, by = points[bi].tolist()
        context.move_to(x0 + ax*scale, y0 + ay*scale)
        context.line_to(x0 + bx*scale, y0 + by*scale)
        context.stroke()
    context.restore()
    # Draw vertices as translucent circles.
    context.save()
    context.set_source_rgba(0.2, 0.2, 1.0, 0.5)
    for point in points:
        x, y = point.tolist()
        nx = x0 + x*scale
        ny = y0 + y*scale
        dot_radius = 2.0
        context.arc(nx, ny, dot_radius, 0, 2*math.pi)
        context.fill()
    context.restore()
    # create the image
    return cairo_helper.get_image_string()

def D_to_DX_internal_first(D, q):
    """
    In the distance matrix the internal vertices come before the leaves
    @param D: the distance matrix with ordered rows and columns
    @param q: the number of leaves
    """
    p = len(D) - q
    return D[:p, p:]

def D_to_DQ_internal_first(D, q):
    """
    In the distance matrix the leaves come before the internal vertices
    @param D: the distance matrix with ordered rows and columns
    @param q: the number of leaves
    """
    p = len(D) - q
    return D[p:, p:]

def D_to_DX_leaves_first(D, q):
    """
    In the distance matrix the leaves come before the internal vertices
    @param D: the distance matrix with ordered rows and columns
    @param q: the number of leaves
    """
    return D[:q, q:]

def D_to_DU_leaves_first(D, q):
    """
    In the distance matrix the leaves come before the internal vertices
    @param D: the distance matrix with ordered rows and columns
    @param q: the number of leaves
    """
    return D[q:, q:]

def D_to_D_star(D, q):
    return D[:q, :q]

def get_column_mean_matrix(M):
    nrows = len(M)
    row = np.mean(M, axis=0)
    R = np.vstack([row]*M)
    if R.shape != M.shape:
        raise ValueError('internal shape fail')
    return R

def get_row_mean_matrix(M):
    R = get_column_mean_matrix(M.T).T
    if R.shape != M.shape:
        raise ValueError('internal shape fail')
    return R

def D_to_B_third(D, q):
    """
    This is the third attempt.
    It assumes that internal nodes are first.
    """
    p = len(D) - q
    DX = D_to_DX_internal_first(D, q)
    DQ = D_to_DQ_internal_first(D, q)
    DX_row_mean_matrix = get_row_mean_matrix(DX)
    DQ_column_mean_matrix = np.vstack([get_column_mean_matrix(DQ)[0]]*p)
    DQ_grand_mean = np.mean(DQ) * np.ones_like(DX)
    try:
        B = DX - DX_row_mean_matrix - DQ_column_mean_matrix + DQ_grand_mean
    except ValueError, e:
        arr = [
                DX.shape,
                DX_row_mean_matrix.shape,
                DQ_column_mean_matrix.shape,
                DQ_grand_mean.shape]
        msg = ', '.join(str(x) for x in arr)
        raise ValueError(msg)
    return B


def D_to_B(D, q):
    """
    This was the first attempt.
    It assumed that leaves were first.
    """
    DX = D_to_DX(D, q)
    DU = D_to_DU(D, q)
    DX_row_mean_matrix = get_row_mean_matrix(DX)
    DU_column_mean_matrix = np.vstack([get_column_mean_matrix(DU)[0]]*q)
    DU_grand_mean = np.mean(DU) * np.ones_like(DX)
    try:
        B = DX - DX_row_mean_matrix - DU_column_mean_matrix - DU_grand_mean
    except ValueError, e:
        arr = [
                DX.shape,
                DX_row_mean_matrix.shape,
                DU_column_mean_matrix.shape,
                DU_grand_mean.shape]
        msg = ', '.join(str(x) for x in arr)
        raise ValueError(msg)
    return B

def D_to_B_second(D, q):
    """
    I think this assumed that leaves were first
    """
    p = len(D)-q
    DX = D_to_DX(D, q)
    DS = D_to_D_star(D, q)
    row_mean_DS = np.mean(DS, axis=1)
    row_mean_matrix = np.vstack([row_mean_DS]*p).T
    col_mean_matrix = get_column_mean_matrix(DX)
    try:
        B = DX - col_mean_matrix - row_mean_matrix + np.mean(DS)
    except ValueError, e:
        arr = [
                DX.shape,
                col_mean_matrix.shape,
                row_mean_matrix.shape]
        msg = ', '.join(str(x) for x in arr)
        raise ValueError(msg)
    return B

if __name__ == '__main__':
    unittest.main()
