"""
Compare double centering vs hollowing of a positive semidefinite matrix.
"""

from StringIO import StringIO
import argparse

import numpy as np
import scipy
from scipy import linalg

import Form
import FormOut
import combobreaker

def get_form():
    """
    @return: the body of a form
    """
    # define the form objects
    form_objects = [
            Form.Integer('length', 'vector length', 5, low=2, high=10)]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_centering_vector(M):
    """
    @param M: symmetric matrix
    """
    return 0.5*np.mean(M) - np.mean(M, axis=0)

def get_hollowing_vector(M):
    """
    @param M: symmetric matrix
    """
    return -0.5*np.diag(M)

def get_delta(v):
    """
    @param v: a centering or hollowing vector
    @return: a matrix to be added to perform the centering or hollowing
    """
    e = np.ones_like(v)
    return np.outer(e, v) + np.outer(v, e)

def get_response_content(fs):
    np.set_printoptions(linewidth=200)
    out = StringIO()
    # get a matrix full of random entries and get the gramian psd matrix
    N = fs.length
    X = np.random.randn(N, N)
    M = np.dot(X.T, X)
    # Hollow out S in one step.
    # hv is the hollowing vector
    # hm is (additive) hollowing matrix
    # ac means it is skipping intermediate step b
    hv_ac = get_hollowing_vector(M)
    hm_ac = get_delta(hv_ac)
    M_ac = M + hm_ac
    # Hollow out S in two steps.
    # d means double centering
    dv_ab = get_centering_vector(M)
    dm_ab = get_delta(dv_ab)
    M_ab = M + dm_ab
    hv_bc = get_hollowing_vector(M_ab)
    hm_bc = get_delta(hv_bc)
    M_bc = M_ab + hm_bc
    # print the results
    print >> out, 'psd (gramian) matrix M:'
    print >> out, M
    print >> out
    print >> out, 'eigenvalues of M:'
    print >> out, scipy.linalg.eigvalsh(M)
    print >> out
    print >> out, 'single step hollowing vector:'
    print >> out, hv_ac
    print >> out
    print >> out, 'single step additive hollowing matrix:'
    print >> out, hm_ac
    print >> out
    print >> out, 'M hollowed out in a single step:'
    print >> out, M_ac
    print >> out
    print >> out, 'first step double centering vector:'
    print >> out, dv_ab
    print >> out
    print >> out, 'first step additive double centering matrix:'
    print >> out, dm_ab
    print >> out
    print >> out, 'double centered M:'
    print >> out, M_ab
    print >> out
    print >> out, 'second step hollowing vector:'
    print >> out, hv_bc
    print >> out
    print >> out, 'second step hollowing matrix:'
    print >> out, hm_bc
    print >> out
    print >> out, 'M hollowed out in two steps:'
    print >> out, M_bc
    print >> out
    return out.getvalue()

