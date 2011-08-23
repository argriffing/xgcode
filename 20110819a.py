"""
Check an identity involving Schur complements, submatrices, and pseudoinverse.
"""

from StringIO import StringIO
import numpy as np
import scipy
import scipy.linalg

import Form
import FormOut

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
            Form.RadioGroup('options', 'symmetry of M', [
                Form.RadioItem('symmetric', 'symmetric', True),
                Form.RadioItem('asymmetric', 'asymmetric')])]
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
    return np.dot(H, np.linalg.inv(np.dot(M, H) + P))

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
    HMH = np.dot(H, np.dot(M, H))
    return np.linalg.inv(HMH + P) - P

def schur(M, nsmall):
    B = M[:nsmall, nsmall:]
    C = np.linalg.inv(M[nsmall:, nsmall:])
    return M[:nsmall, :nsmall] - np.dot(B, np.dot(C, B.T))

def assert_named_equation(a_name_pair, b_name_pair):
    """
    This is a helper function to check for counterexamples.
    """
    a, a_name = a_name_pair
    b, b_name = b_name_pair
    if not np.allclose(a, b):
        out = StringIO()
        print >> out, a_name + ':'
        print >> out, a
        print >> out, b_name + ':'
        print >> out, b
        raise Counterexample(out.getvalue())

def assert_pinvproj(fs, M):
    """
    Raise a Counterexample if one is found.
    """
    pinvproj_of_sub = pinvproj(M[:fs.nsmall, :fs.nsmall])
    schur_of_pinvproj = schur(pinvproj(M), fs.nsmall)
    bottduff_of_sub = bott_duffin(M[:fs.nsmall, :fs.nsmall])
    schur_of_bottduff = schur(bott_duffin(M), fs.nsmall)
    assert_named_equation(
            (pinvproj_of_sub, 'pinvproj of sub'),
            (schur_of_pinvproj, 'schur of pinvproj'))
    assert_named_equation(
            (pinvproj_of_sub, 'pinvproj of sub'),
            (bottduff_of_sub, 'bottduff of sub'))
    assert_named_equation(
            (schur_of_pinvproj, 'schur of pinvproj'),
            (schur_of_bottduff, 'schur of bottduff'))

def double_centered(M):
    """
    This is slow.
    """
    nrows, ncols = M.shape
    if nrows != ncols:
        raise ValueError('expected a square matrix')
    e = np.ones(nrows)
    I = np.eye(nrows)
    P = np.outer(e, e) / np.inner(e, e)
    H = I - P
    return np.dot(H, np.dot(M, H))

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
    HMH = double_centered(M)
    HMH_A = HMH[:fs.nsmall, :fs.nsmall]
    # check an identity between two ways to compute a centered submatrix
    HAH_direct = double_centered(A)
    HAH_indirect = double_centered(HMH_A)
    assert_named_equation(
            (HAH_direct, 'double centered submatrix'),
            (HAH_indirect, 're-centered submatrix of doubly centered matrix'))
    HAH = HAH_indirect
    # This is not true:
    # check that HAH <==inverse_in_H==> H A^-1 H
    HAH_pinv = inverse_in_H(HAH)
    H_Ainv_H = double_centered(np.linalg.inv(A))
    assert_named_equation(
            (HAH_pinv, 'inverse-in-H of HAH'),
            (H_Ainv_H, 'double centered inverse of A'))
    # check a submatrix-schur-inverse commutativity in double centering
    HMH_pinv = inverse_in_H(HMH)
    schur_of_HMH_pinv = schur(HMH_pinv, fs.nsmall)
    assert_named_equation(
            (HAH_pinv, 'inverse-in-H of HAH'),
            (schur_of_HMH_pinv, 'schur of pinv of HMH'))



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

def MatrixSampler:
    def __init__(self, nbig):
        self.nbig = nbig

def SymmetricMatrixSampler(MatrixSampler):
    pass

def AsymmetricMatrixSampler(MatrixSampler):
    # create a random matrix
    B = 10 * (np.random.rand(fs.nbig, fs.nbig) - 0.5)
    if fs.symmetric:
        M = B + B.T
    else:
        M = B

def DistanceMatrixSampler(MatrixSampler):
    def __call__(self):
        """
        There is probably a nicer way to make an EDM.
        """
        nrows = self.nbig
        ncols = self.nbig-1
        X = np.random.rand(nrows, ncols)
        M = np.zeros((self.nbig, self.nbig))
        for i in range(self.nbig):
            for j in range(self.nbig):
                d = X[i] - X[j]
                M[i, j] = np.dot(d, d)
        return M

def search_for_counterexample(fs, sampler):
    """
    Raise a Counterexample if one is found.
    """
    for i in range(nsamples):
        # draw a matrix from the distribution
        M = sampler()
        # assert pinvproj and bott-duffin identities
        assert_pinvproj(fs, M)
        # assert identities regarding double centering
        assert_double_centering_identities(fs, M)


def search_for_generic_counterexample(fs):
    for i in range(fs.nsamples):
        # create a random matrix
        B = 10 * (np.random.rand(fs.nbig, fs.nbig) - 0.5)
        if fs.symmetric:
            M = B + B.T
        else:
            M = B
        # assert pinvproj and bott-duffin identities
        assert_pinvproj(fs, M)
        # assert identities regarding double centering
        assert_double_centering_identities(fs, M)

def search_for_distance_counterexample(fs):
    """
    Raise a Counterexample if one is found.
    """
    for i in range(fs.nsamples):
        # create a random matrix
        M = sample_edm(fs.nbig)
        # assert pinvproj and bott-duffin identities
        assert_pinvproj(fs, M)
        # assert identities regarding double centering
        assert_double_centering_identities(fs, M)

def get_response_content(fs):
    out = StringIO()
    try:
        search_for_distance_counterexample(fs)
    except Counterexample, e:
        print >> out, 'Found a counterexample to a conjectured identity'
        print >> out, 'using a random distance matrix.'
        print >> out, str(e)
        return out.getvalue()
    print >> out, 'no counterexample found using random distance matrices'
    try:
        search_for_generic_counterexample(fs)
    except Counterexample, e:
        print >> out, 'Found a counterexample to a conjectured identity'
        print >> out, 'using a random generic matrix.'
        print >> out, str(e)
        return out.getvalue()
    print >> out, 'no counterexample found using random generic matrices'
    return out.getvalue()

