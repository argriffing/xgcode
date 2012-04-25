"""
Matrix formatting tools.

This module uses numpy arrays and not numpy matrices.
"""

import unittest
from StringIO import StringIO

import numpy as np

import matrixio
import graph


#TODO put these global arrays into const data files

# this is a rotation of a contrast matrix of the example tree in "Why neighbor joining works"
g_example_loading_matrix = np.array([
    [0.187795256419, -0.178119343731, 0.107859503348, -0.566026656126, -0.561620966297, 1.71570965693, -2.65637730248],
    [0.187795256419, -0.178119343731, 0.107859503348, -0.566026656126, -0.561620966297, -1.71570965693, 2.65637730248],
    [0.181174971619, -0.156327719297, 0.0923993967909, 1.1689340078, 1.17504869278, 0.0, 0.0],
    [0.0942074702731, 0.512566406759, -0.607716820144, -0.0368806955525, -0.0206980063065, 0.0, 0.0],
    [-0.187795256419, -0.178119343731, -0.107859503348, -0.566026656126, 0.561620966296, 2.65637730248, 1.71570965693],
    [-0.187795256419, -0.178119343731, -0.107859503348, -0.566026656126, 0.561620966296, -2.65637730248, -1.71570965693],
    [-0.181174971619, -0.156327719297, -0.0923993967909, 1.16893400781, -1.17504869278, 0.0, 0.0],
    [-0.0942074702731, 0.512566406759, 0.607716820144, -0.0368806955526, 0.0206980063064, 0.0, 0.0]])

# this is the varimax rotation of the example loading matrix for eps=1e-5
g_example_rotated_matrix = np.array([
    [0.330815132017, -0.256910070995, 0.0801815618365, -0.534892400091, -0.496438433195, -0.000808741457446, -3.16227755675],
    [0.330815132017, -0.256910070995, 0.0801815618365, -0.534892400091, -0.496438433195, 0.000808741457446, 3.16227755675],
    [-0.0786945562714, 0.0117297718391, 0.00988669420354, 1.17928258946, 1.18987853431, 0.0, 0.0],
    [-0.0371881711244, 0.502090370696, -0.612105936609, -0.109497791434, -0.0506611396544, 0.0, 0.0],
    [-0.330815129488, -0.256910072776, -0.0801815614454, -0.534892421301, 0.496438411167, 3.16227755675, -0.000808741457446],
    [-0.330815129488, -0.256910072776, -0.0801815614454, -0.534892421301, 0.496438411167, -3.16227755675, 0.000808741457446],
    [0.0786945522922, 0.0117297739385, -0.00988669475083, 1.17928263783, -1.18987848662, 0.0, 0.0],
    [0.0371881700454, 0.50209037107, 0.612105936374, -0.109497793077, 0.0506611360203, 0.0, 0.0]])

class MatrixError(ValueError): pass

def ndot(*args):
    M = args[0]
    for B in args[1:]:
        M = np.dot(M, B)
    return M

def assert_1d(M):
    if len(M.shape) != 1:
        raise MatrixError('the array is not 1d')

def assert_2d(M):
    if len(M.shape) != 2:
        raise MatrixError('the array is not 2d')

def assert_square(M):
    assert_2d(M)
    nrows, ncols = M.shape
    if nrows != ncols:
        raise MatrixError('the matrix is not square')

def assert_sign_symmetric(M):
    assert_square(M)
    S = np.sign(M).astype(int)
    if not np.array_equal(S, S.T):
        raise MatrixError('the matrix is not sign symmetric')

def assert_symmetric(M):
    assert_square(M)
    if not np.allclose(M, M.T):
        raise MatrixError('the matrix is not symmetric')

def assert_nonnegative(M):
    if (M < 0).any():
        raise MatrixError('the matrix has a negative element')

def assert_hollow(M):
    assert_square(M)
    if np.diag(M).any():
        raise MatrixError('the matrix is not hollow')

def assert_predistance(M):
    assert_symmetric(M)
    assert_nonnegative(M)

def assert_weighted_adjacency(M):
    assert_symmetric(M)
    assert_nonnegative(M)

def assert_symmetric_irreducible(M):
    assert_symmetric(M)
    V = set(range(len(M)))
    E = set()
    for i in range(len(M)):
        for j in range(i):
            if M[i, j]:
                edge = frozenset((i,j))
                E.add(edge)
    nd = graph.g_to_nd(V, E)
    V_component, D_component = graph.nd_to_dag_component(nd, 0)
    if V_component != V:
        raise MatrixError('the matrix is not irreducible')

#FIXME do something about this
def assert_detailed_balance(M):
    pass

def assert_rate_matrix(M):
    """
    @param M: a numpy array
    """
    assert_square(M)
    n = M.shape[0]
    if not (np.diag(M) < 0).all():
        raise MatrixError('diagonal elements of a rate matrix should be less than zero')
    if ((M - np.diag(M)) < 0).any():
        raise MatrixError('off diagonal elements of a rate matrix should be nonnegative')
    if not np.allclose(np.sum(M, 1), np.zeros(n)):
        raise MatrixError('rows of a rate matrix should each sum to zero')

def assert_transition_matrix(M):
    """
    @param M: a numpy array
    """
    assert_square(M)
    n = M.shape[0]
    if (M < 0).any():
        raise MatrixError('each element of a transition matrix should be nonnegative')
    if not np.allclose(np.sum(M, 1), np.ones(n)):
        raise MatrixError('rows of a transition matrix should each sum to 1.0')

def row_major_to_dict(matrix, ordered_row_labels, ordered_column_labels):
    """
    Convert between matrix formats.
    @param matrix: a matrix represented a row major list of lists
    @param ordered_row_labels: the row labels
    @param ordered_column_labels: the column labels
    @return: a matrix represented as a dictionary of (row, column) label pairs
    """
    # do some sanity checking
    nrows = len(ordered_row_labels)
    ncols = len(ordered_column_labels)
    assert matrix
    assert len(matrix) == nrows
    assert len(matrix[0]) == ncols
    # create the matrix
    dict_matrix = {}
    for row_label, row in zip(ordered_row_labels, matrix):
        for column_label, element in zip(ordered_column_labels, row):
            key = (row_label, column_label)
            dict_matrix[key] = element
    return dict_matrix

def dict_to_row_major(matrix, ordered_row_labels, ordered_column_labels):
    """
    Convert between matrix formats.
    @param matrix: a matrix represented as a dictionary of (row, column) label pairs
    @param ordered_row_labels: the row labels
    @param ordered_column_labels: the column labels
    @return: a row major matrix without label information
    """
    # do some sanity checking
    nrows = len(ordered_row_labels)
    ncols = len(ordered_column_labels)
    assert len(matrix) == nrows * ncols
    # create the matrix
    row_major_matrix = []
    for a in ordered_row_labels:
        row = []
        for b in ordered_column_labels:
            row.append(matrix[(a, b)])
        row_major_matrix.append(row)
    return row_major_matrix

def m_to_matlab_string(matrix):
    """
    @param matrix: a list of lists
    @return: a single line string suitable for copying and pasting into matlab
    """
    return '[' + '; '.join(' '.join(str(x) for x in row) for row in matrix) + ']'

def list_to_matrix(arr, f):
    """
    Apply a function to pairwise elements of a list.
    @param arr: a list of values, pairs of which will be passed to f
    @param f: a function that takes two arguments
    @return: a symmetric row major matrix
    """
    #FIXME this function is used but is stupidly named
    M = []
    for a in arr:
        row = [f(a, b) for b in arr]
        M.append(row)
    return M

def get_principal_submatrix(M, indices):
    """
    @param M: a square numpy array
    @param indices: indices of M to be included in the output matrix
    @return: a square numpy array no bigger than M
    """
    assert_square(M)
    n_small = len(indices)
    R = np.zeros((n_small, n_small))
    for i_small, i_big in enumerate(indices):
        for j_small, j_big in enumerate(indices):
            R[i_small][j_small] = M[i_big][j_big]
    return R

def double_centered(M):
    """
    @param M: a numpy array
    @return: a doubly centered numpy array
    """
    assert_square(M)
    n = len(M)
    r = np.mean(M, 0)
    c = np.mean(M, 1)
    m = np.mean(M)
    e = np.ones(n)
    HMH = M - np.outer(e, r) - np.outer(c, e) + m
    return HMH

def varimax(M, eps=1e-5, itermax=200):
    """
    @param M: a numpy matrix to rotate
    @param eps: use smaller values to give a more accurate rotation
    @param itermax: do at most this many iterations
    @return: a rotation matrix R such that MR rotates M to satisfy the varimax criterion
    """
    nrows, ncols = M.shape
    R = np.eye(ncols)
    d = 0
    for i in range(itermax):
        MR = np.dot(M, R)
        D = np.diag(np.sum(MR**2, 0))/nrows
        U, S, V_T = np.linalg.svd(np.dot(M.T, MR**3 - np.dot(MR, D)))
        R = np.dot(U, V_T)
        d_next = np.sum(S)
        if d_next < d * (1 + eps):
            break
        d = d_next
    return R

def get_best_reflection(A, B):
    """
    Find best {-1, +1} multiples of columns of A to approximate B.
    That is, minimize the frobenius norm of AS - B
    where S is a diagonal matrix with {-1, +1} diagonal entries.
    The returned array has the diagonal elements of S.
    A motivation for this function is for the rows of matrix A
    to represent new points up to reflection across each axis,
    whereas rows of matrix B represent old points.
    If v is the output of this function then, using numpy notation,
    A*v gives the canonically reflected points of A
    that best map onto the points of B.
    @param A: a matrix
    @param B: a matrix
    @return: a 1d array of {-1, +1} entries
    """
    assert_2d(A)
    assert_2d(B)
    if A.shape != B.shape:
        raise MatrixError('the shapes of A and B should be the same')
    nrows, ncols = A.shape
    best_signs = np.ones(ncols)
    for i in range(ncols):
        pos_error = np.linalg.norm(A.T[i] - B.T[i])
        neg_error = np.linalg.norm(-A.T[i] - B.T[i])
        if neg_error < pos_error:
            best_signs[i] = -1
    return best_signs

class TestMatrixUtil(unittest.TestCase):

    def test_sign_symmetric_true(self):
        M = np.array([
                [0.0, -0.5],
                [-0.5, 0.0]])
        assert_sign_symmetric(M)

    def test_sign_symmetric_false(self):
        M = np.array([
                [0.0, 0.5],
                [-0.5, 0.0]])
        self.assertRaises(MatrixError, assert_sign_symmetric, M)

    def test_symmetric_irreducible_true(self):
        M = np.array([
                [1, 1, 0],
                [1, 1, 1],
                [0, 1, 1]])
        assert_symmetric_irreducible(M)

    def test_symmetric_irreducible_false_reducible(self):
        M = np.array([
                [1, 0, 0],
                [0, 1, 1],
                [0, 1, 1]])
        self.assertRaises(MatrixError, assert_symmetric_irreducible, M)

    def test_symmetric_irreducible_false_asymmetric(self):
        M = np.array([
                [0, 1],
                [0, 0]])
        self.assertRaises(MatrixError, assert_symmetric_irreducible, M)

    def test_dict_to_row_major(self):
        d = {
                (0, 0) : 1, (0, 1) : 2, (0, 2) : 3,
                (1, 0) : 4, (1, 1) : 5, (1, 2) : 6
                }
        expected = [[1, 2, 3], [4, 5, 6]]
        row_labels = [0, 1]
        column_labels = [0, 1, 2]
        observed = dict_to_row_major(d, row_labels, column_labels)
        self.assertEquals(expected, observed)

    def test_read_square_matrix_bad(self):
        s = '1 2 3 \n 4 5 6'
        M = np.array(matrixio.read_matrix(StringIO(s)))
        self.assertRaises(MatrixError, assert_square, M)

    def test_read_square_matrix_good(self):
        s = '1 2 3 \n 4 5 6 \n 7 8 9'
        row_major_observed = matrixio.read_matrix(StringIO(s))
        row_major_expected = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
        self.assertEquals(row_major_observed, row_major_expected)

    def test_read_rate_matrix(self):
        arr = [[-3, 1, 1, 1], [1, -3, 1, 1], [1, 1, -3, 1], [1, 1, 1, -3]]
        s = '\n'.join('\t'.join(str(x) for x in row) for row in arr)
        sio = StringIO(s)
        M = np.array(matrixio.read_matrix(sio))
        assert_rate_matrix(M)

    def test_m_to_matlab_string(self):
        observed = m_to_matlab_string([[1, 2, 3], [4, 5, 6]])
        expected = '[1 2 3; 4 5 6]'
        self.assertEquals(expected, observed)

    def test_varimax(self):
        M = g_example_loading_matrix
        MR = np.dot(M, varimax(M))
        self.assertTrue(np.allclose(MR, g_example_rotated_matrix))

    def test_double_centered_a(self):
        M = np.array([[1.0, 2.0], [3.0, 4.0]])
        observed = double_centered(M)
        expected = np.array([[0, 0], [0, 0]])
        self.assertTrue(np.allclose(observed, expected))

    def test_double_centered_b(self):
        M = np.array([[1.0, 2.0], [1.0, 1.0]])
        observed = double_centered(M)
        expected = np.array([[-0.25, 0.25], [0.25, -0.25]])
        self.assertTrue(np.allclose(observed, expected))


if __name__ == '__main__':
    unittest.main()

