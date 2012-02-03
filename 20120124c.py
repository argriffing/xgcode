"""
Does adding selection always decrease the expected rate.
"""

from StringIO import StringIO
import math

import numpy as np
import scipy
from scipy import linalg

import Form
import FormOut
from MatrixUtil import ndot
import mrate
import divtime
import ctmcmi

def get_form():
    """
    @return: the body of a form
    """
    # define the form objects
    form_objects = [
            Form.Float('p_mid',
                'mutation process probability for intermediate state',
                '0.01', low_exclusive=0, high_exclusive=1)]
            #Form.Integer('nstates', 'number of states', 4, low=2, high=9)]
            #Form.Float('divtime', 'arbitrary large-ish divergence time',
                #'3', low_exclusive=0)]
            #Form.Float('delta', 'vanishing time delta',
                #'0.0001', high_exclusive=0.25)]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_barbell_rate_matrix(p_mid):
    # define a hollow exchangeability-like matrix
    nstates = 3
    Z = np.array([
        [0.0, 1.0, 0.0],
        [1.0, 0.0, 1.0],
        [0.0, 1.0, 0.0]])
    # define the stationary distribution
    p = np.array([(1 - p_mid)/2, p_mid, (1 - p_mid)/2])
    # define the mutation matrix
    D = np.diag(p)
    D_inv = np.diag(np.reciprocal(p))
    Q_unnormal = ndot(D_inv**0.5, Z, D**0.5)
    Q = np.copy(Q_unnormal)
    for i in range(nstates):
        Q[i, i] = -np.sum(Q[i])
    return Q, p

def get_response_content(fs):
    out = StringIO()
    np.set_printoptions(linewidth=200)
    # define the barbell mutation rate matrix
    M, p = get_barbell_rate_matrix(fs.p_mid)
    nstates = len(p)
    print >> out, 'barbell mutation matrix:'
    print >> out, M
    print >> out
    print >> out, 'all of these should be zero for detailed balance:'
    for i in range(nstates):
        for j in range(nstates):
            print >> out, p[i] * M[i, j] - p[j]*M[j, i]
    print >> out
    print >> out, 'expected rate of the barbell mutation matrix:'
    print >> out, mrate.Q_to_expected_rate(M)
    print >> out
    p_target = np.array([1/3., 1/3., 1/3.])
    print >> out, 'target stationary distribution:'
    print >> out, p_target
    print >> out
    Q = divtime.to_gtr_halpern_bruno(M, p_target)
    print >> out, 'mutation-selection balance rate matrix:'
    print >> out, Q
    print >> out
    v = mrate.R_to_distn(Q)
    print >> out, 'computed stationary distribution:'
    print >> out, v
    print >> out
    print >> out, 'expected rate of the mutation-selection balance rate matrix:'
    print >> out, mrate.Q_to_expected_rate(Q)
    print >> out
    print >> out, 'all of these should be zero for detailed balance:'
    for i in range(nstates):
        for j in range(nstates):
            print >> out, v[i] * Q[i, j] - v[j]*Q[j, i]
    print >> out
    return out.getvalue()

