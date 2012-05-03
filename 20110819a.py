"""
Check an identity involving Schur complements, submatrices, and pseudoinverse.
"""

from StringIO import StringIO

import numpy as np
import scipy
import scipy.linalg

import Form
import FormOut
import numpyutils
import MatrixUtil
from MatrixUtil import double_centered_slow
from MatrixUtil import ndot

# TODO use combobreaker for the counterexample search

class Counterexample(Exception): pass

def get_form():
    """
    @return: the body of a form
    """
    # define the form objects
    form_objects = [
            Form.Integer('nsamples', 'number of samples',
                1000, low=1, high=10000),
            Form.Integer('nbig', 'order of M',
                10, low=1, high=12),
            Form.Integer('nsmall', 'order of principal submatrix',
                4, low=1, high=12),
            Form.RadioGroup('matrix_class', 'matrix class', [
                Form.RadioItem('distance', 'distance', True),
                Form.RadioItem('predistance', 'predistance'),
                Form.RadioItem('nonnegsym', 'nonnegative symmetric'),
                Form.RadioItem('symmetric', 'symmetric'),
                Form.RadioItem('asymmetric', 'asymmetric'),
                Form.RadioItem('posdef', 'positive definite'),
                Form.RadioItem('negdef', 'negative definite')])]
    return form_objects

def get_form_out():
    return FormOut.Report()

def bott_duffin(M):
    """
    We pretend that P_L is H.
    """
    nrows, ncols = M.shape
    if nrows != ncols:
        raise ValueError('expected a square matrix')
    e = np.ones(nrows)
    I = np.eye(nrows)
    P = np.outer(e, e) / np.inner(e, e)
    H = I - P
    return ndot(H, np.linalg.inv(np.dot(M, H) + P))

def pinvproj(M):
    """
    This could be made more efficient using double centering.
    """
    nrows, ncols = M.shape
    if nrows != ncols:
        raise ValueError('expected a square matrix')
    e = np.ones(nrows)
    I = np.eye(nrows)
    P = np.outer(e, e) / np.inner(e, e)
    H = I - P
    #
    #HMH = MatrixUtil.double_centered(M)
    #return np.linalg.pinv(HMH)
    #
    HMH = ndot(H, M, H)
    W = HMH + P
    d = np.linalg.det(W)
    if abs(d) < 1e-5:
        #raise ValueError('small determinant in pinvproj: %f' % d)
        pass
    return np.linalg.inv(W) - P

def schur(M, nsmall):
    A = M[:nsmall, :nsmall]
    B = M[:nsmall, nsmall:]
    C = M[nsmall:, nsmall:]
    d = np.linalg.det(C)
    if abs(d) < 1e-5:
        raise ValueError('small determinant for schur complement')
    C_inv = np.linalg.inv(C)
    return A - ndot(B, C_inv, B.T)

def assert_named_equation(a_name_pair, b_name_pair, M):
    """
    This is a helper function to check for counterexamples.
    Check that matrix A is the same as matrix B.
    If they are not the same, report M and its eigendecomposition.
    @param a_name_pair: (matrix A, name of A)
    @param b_name_pair: (matrix B, name of B)
    @param M: underlying matrix M
    """
    a, a_name = a_name_pair
    b, b_name = b_name_pair
    if not np.allclose(a, b):
        w, vt = scipy.linalg.eigh(M)
        out = StringIO()
        print >> out, a_name + ':'
        print >> out, a
        print >> out, b_name + ':'
        print >> out, b
        print >> out, 'M:'
        print >> out, M
        print >> out, 'eigenvalues of M:'
        print >> out, w
        print >> out, 'columns are eigenvectors of M:'
        print >> out, vt
        raise Counterexample(out.getvalue())

def assert_pinvproj(fs, M):
    """
    Raise a Counterexample if one is found.
    """
    M_sub = M[:fs.nsmall, :fs.nsmall]
    pinvproj_of_sub = pinvproj(M_sub)
    schur_of_pinvproj = schur(pinvproj(M), fs.nsmall)
    bottduff_of_sub = numpyutils.bott_duffin_const(M_sub)
    schur_of_bottduff = schur(numpyutils.bott_duffin_const(M), fs.nsmall)
    """
    assert_named_equation(
            (np.array([[1]]), 'one'),
            (np.array([[2]]), 'two'),
            MatrixUtil.double_centered(M))
    """
    assert_named_equation(
            (pinvproj_of_sub, 'pinvproj of sub'),
            (schur_of_pinvproj, 'schur of pinvproj'),
            double_centered_slow(M))
    assert_named_equation(
            (pinvproj_of_sub, 'pinvproj of sub'),
            (bottduff_of_sub, 'bottduff of sub'),
            double_centered_slow(M))
    assert_named_equation(
            (schur_of_pinvproj, 'schur of pinvproj'),
            (schur_of_bottduff, 'schur of bottduff'),
            double_centered_slow(M))

def inverse_in_H(M):
    nrows, ncols = M.shape
    if nrows != ncols:
        raise ValueError('expected a square matrix')
    e = np.ones(nrows)
    I = np.eye(nrows)
    P = np.outer(e, e) / np.inner(e, e)
    return np.linalg.inv(M + P) - P

def assert_double_centering_identities(fs, M):
    """
    Raise a Counterexample if one is found.
    """
    A = M[:fs.nsmall, :fs.nsmall]
    HMH = double_centered_slow(M)
    HMH_A = HMH[:fs.nsmall, :fs.nsmall]
    # check an identity between two ways to compute a centered submatrix
    HAH_direct = double_centered_slow(A)
    HAH_indirect = double_centered_slow(HMH_A)
    assert_named_equation(
            (HAH_direct, 'double centered submatrix'),
            (HAH_indirect, 're-centered submatrix of doubly centered matrix'),
            M)
    HAH = HAH_indirect
    # This is not true:
    # check that HAH <==inverse_in_H==> H A^-1 H
    HAH_pinv = inverse_in_H(HAH)
    H_Ainv_H = double_centered_slow(np.linalg.inv(A))
    """
    assert_named_equation(
            (HAH_pinv, 'inverse-in-H of HAH'),
            (H_Ainv_H, 'double centered inverse of A'))
    """
    # check a submatrix-schur-inverse commutativity in double centering
    HMH_pinv = inverse_in_H(HMH)
    schur_of_HMH_pinv = schur(HMH_pinv, fs.nsmall)
    assert_named_equation(
            (HAH_pinv, 'inverse-in-H of HAH'),
            (schur_of_HMH_pinv, 'schur of pinv of HMH'),
            M)



def wat():
    print >> out, 'nsamples:', fs.nsamples
    print >> out, 'no counterexamples found'
    print >> out, 'example identity:'
    print >> out, 'M:'
    print >> out, M
    print >> out, 'pinvproj of sub-matrix:'
    print >> out, pinvproj_of_sub
    print >> out, 'schur complement of pinvproj:'
    print >> out, schur_of_pinvproj
    print >> out, 'bottduff of sub-matrix:'
    print >> out, bottduff_of_sub
    print >> out, 'schur complement of bottduff:'
    print >> out, schur_of_bottduff

def sample_square_matrix(n):
    """
    Sample a very generic square matrix.
    """
    return 10.0 * np.random.rand(n, n) - 5.0

class MatrixSampler:
    def __init__(self, n):
        self.n = n

class DMSampler(MatrixSampler):
    def __call__(self):
        nrows = self.n
        ncols = self.n-1
        X = np.random.rand(nrows, ncols)
        M = np.zeros((self.n, self.n))
        for i in range(self.n):
            for j in range(self.n):
                d = X[i] - X[j]
                M[i, j] = np.dot(d, d)
        return M

class PredistanceMatrixSampler(MatrixSampler):
    def __call__(self):
        B = sample_square_matrix(self.n)
        S = np.abs(B + B.T)
        return S - np.diag(np.diag(S))

class NonNegSymmetricMatrixSampler(MatrixSampler):
    def __call__(self):
        B = sample_square_matrix(self.n)
        S = np.abs(B + B.T)
        return S

class SymmetricMatrixSampler(MatrixSampler):
    def __call__(self):
        B = sample_square_matrix(self.n)
        return B + B.T

class AsymmetricMatrixSampler(MatrixSampler):
    def __call__(self):
        return sample_square_matrix(self.n)

class PositiveDefiniteMatrixSampler(MatrixSampler):
    def __call__(self):
        B = sample_square_matrix(self.n)
        return np.dot(B, B.T)

class NegativeDefiniteMatrixSampler(MatrixSampler):
    def __call__(self):
        B = sample_square_matrix(self.n)
        return -np.dot(B, B.T)

def search_for_counterexample(fs, sampler):
    """
    Raise a Counterexample if one is found.
    """
    for i in range(fs.nsamples):
        # draw a matrix from the distribution
        M = sampler()
        # assert pinvproj and bott-duffin identities
        assert_pinvproj(fs, M)
        # assert identities regarding double centering
        assert_double_centering_identities(fs, M)

def get_response_content(fs):
    # tell numpy to print verbosely
    np.set_printoptions(linewidth=200)
    # select the matrix sampler chosen by the user
    if fs.distance:
        sampler = DMSampler(fs.nbig)
    elif fs.predistance:
        sampler = PredistanceMatrixSampler(fs.nbig)
    elif fs.nonnegsym:
        sampler = NonNegSymmetricMatrixSampler(fs.nbig)
    elif fs.symmetric:
        sampler = SymmetricMatrixSampler(fs.nbig)
    elif fs.asymmetric:
        sampler = AsymmetricMatrixSampler(fs.nbig)
    elif fs.posdef:
        sampler = PositiveDefiniteMatrixSampler(fs.nbig)
    elif fs.negdef:
        sampler = NegativeDefiniteMatrixSampler(fs.nbig)
    else:
        raise ValueError('invalid matrix class')
    # look for a counterexample
    out = StringIO()
    try:
        search_for_counterexample(fs, sampler)
    except Counterexample, e:
        print >> out, 'Found a counterexample:'
        print >> out, str(e)
    else:
        print >> out, 'no counterexample found'
    return out.getvalue()

