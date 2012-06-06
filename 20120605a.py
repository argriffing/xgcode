"""
Check the spectral effect of adding a sink state to a reversible Markov model.
"""

from StringIO import StringIO
import random
import math
import itertools
from itertools import combinations
from itertools import product

import numpy as np
import scipy
from scipy import linalg

import Form
import FormOut
import mrate
import ctmcmi
import cheeger
import msimpl
import iterutils
import MatrixUtil
from MatrixUtil import ndot


def get_form():
    form_objects = [
            Form.Integer('nstates', 'number of states', 4, low=2, high=10),
            ]
    return form_objects

def get_form_out():
    return FormOut.Report()

def process(fs):
    n = fs.nstates
    np.set_printoptions(linewidth=200)
    out = StringIO()
    # Sample a symmetric rate matrix and a stationary distribution,
    # then construct the rate matrix R.
    S = MatrixUtil.sample_pos_sym_matrix(n)
    v = mrate.sample_distn(n)
    psi = np.sqrt(v)
    R = (S.T / psi).T * psi
    R -= np.diag(np.sum(R, axis=1))
    R_W, R_V = scipy.linalg.eig(R)
    # Add some deletion.
    deletion = 0.1
    Q = np.zeros((n+1, n+1))
    for i in range(n):
        for j in range(n):
            Q[i, j] = R[i, j]
    for i in range(n):
        Q[i, -1] = deletion
    Q -= np.diag(np.sum(Q, axis=1))
    Q_W, Q_V = scipy.linalg.eig(Q)
    print >> out, 'deletion rate:'
    print >> out, deletion
    print >> out
    print >> out, 'sampled rate matrix R:'
    print >> out, R
    print >> out
    print >> out, 'spectrum of R:'
    print >> out, R_W
    print >> out
    print >> out, 'rate matrix with deletion Q:'
    print >> out, Q
    print >> out
    print >> out, 'spectrum of Q:'
    print >> out, Q_W
    print >> out
    return out.getvalue().rstrip()

def get_response_content(fs):
    return process(fs) + '\n'

