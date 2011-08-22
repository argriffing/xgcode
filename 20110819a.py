"""
Check an identity involving Schur complements, submatrices, and pseudoinverse.
"""

from StringIO import StringIO
import numpy as np
import scipy
import scipy.linalg
import sympy

import Form
import FormOut

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

def get_response_content(fs):
    out = StringIO()
    for i in range(fs.nsamples):
        # create a random symmetric matrix
        B = 10 * (np.random.rand(fs.nbig, fs.nbig) - 0.5)
        if fs.symmetric:
            M = B + B.T
        else:
            M = B
        pinvproj_of_sub = pinvproj(M[:fs.nsmall, :fs.nsmall])
        schur_of_pinvproj = schur(pinvproj(M), fs.nsmall)
        bottduff_of_sub = bott_duffin(M[:fs.nsmall, :fs.nsmall])
        schur_of_bottduff = schur(bott_duffin(M), fs.nsmall)
        fail = False
        if not np.allclose(pinvproj_of_sub, schur_of_pinvproj):
            print >> out, 'pinvproj counterexample (unexpected!):'
            print >> out, pinvproj_of_sub
            print >> out, schur_of_pinvproj
            fail = True
        if not np.allclose(pinvproj_of_sub, bottduff_of_sub):
            print >> out, 'pinvproj of sub vs. bottduff of sub counter:'
            print >> out, pinvproj_of_sub
            print >> out, bottduff_of_sub
            fail = True
        if not np.allclose(schur_of_pinvproj, schur_of_bottduff):
            print >> out, 'schur of pinvproj vs. schur of bottduff counter:'
            print >> out, schur_of_pinvproj
            print >> out, schur_of_bottduff
            fail = True
        if fail:
            return out.getvalue()
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
    return out.getvalue()

