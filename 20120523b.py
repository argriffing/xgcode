"""
Compute commute times for reversible Markov processes.

For irreducible processes the matrix of commute times is
a spherical EDM in the sense that its NxN entries are
squared Euclidean distances between points in at least N-1
dimensional Euclidean space.
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
    # construct the symmetric matrix that is similar to R
    R_sim = (R.T * psi).T / psi
    if not np.allclose(R_sim, R_sim.T):
        raise ValueError('the similar symmetric matrix is not symmetric...')
    R_sim_W, R_sim_V = scipy.linalg.eigh(R_sim)
    R_gap = -R_sim_W[-2]
    v2 = R_sim_V.T[-2]**2
    # reconstruct the eigenvectors of R
    R_V_rebuilt = (R_sim_V.T / psi).T
    # Try to make the commute time matrix.
    # R_sim is a lot like a Laplacian matrix, so lets pseudoinvert it.
    R_sim_pinv = scipy.linalg.pinv(R_sim)
    myouter = np.outer(np.ones(n), np.diag(R_sim_pinv))
    D = -(myouter + myouter.T - 2*R_sim_pinv)
    D_commute = mrate.get_commute_distance_matrix(R, v)
    if not np.allclose(D, D_commute):
        raise ValueError('error computing commute distances')
    HDH = MatrixUtil.double_centered(D)
    HDH_W, HDH_V = scipy.linalg.eigh(HDH)
    # compute squared pairwise distances brutely
    X = R_sim_V.T[:-1].T / np.sqrt(-R_sim_W[:-1])
    D_brute = np.array([[np.dot(b - a, b - a) for a in X] for b in X])
    print >> out, 'reconstructed EDM:'
    print >> out, D
    print >> out
    D = (D.T / psi).T / psi
    print >> out, 'divide by square roots of stationary probabilities:'
    print >> out, D
    print >> out
    print >> out, 'eigh of centered EDM:'
    print >> out, 'eigenvalues:'
    print >> out, HDH_W
    print >> out, 'reciprocal nonzero eigenvalues:'
    print >> out, 1 / HDH_W
    print >> out, 'eigenvectors:'
    print >> out, HDH_V
    print >> out
    print >> out, 'squared distances computed brutely:'
    print >> out, D_brute
    print >> out
    print >> out, '1 / (h * max(D)):', 1 / (np.dot(v, 1-v) * np.max(D))
    print >> out, '1 / max(D):', 1 / np.max(D)
    print >> out
    # report some more standard stuff
    print >> out, 'sampled rate matrix R:'
    print >> out, R
    print >> out, 'stationary distn:', v
    print >> out, '1/R01 + 1/R10:', 1/R[0,1] + 1/R[1,0]
    print >> out
    print >> out, 'scipy.linagl.eig(R):'
    print >> out, R_W
    print >> out, R_V
    print >> out
    print >> out, 'symmetric matrix similar to R:'
    print >> out, R_sim
    print >> out
    print >> out, 'eigh of the symmetric similar matrix to R:'
    print >> out, R_sim_W
    print >> out, R_sim_V
    print >> out, 'spectral gap:', R_gap
    print >> out, 'entrywise squares of eigenvectors:'
    print >> out, R_sim_V ** 2
    print >> out, 'a bilinear form involving a fiedler-like eigenvector:'
    print >> out, ndot(R_sim_V.T[-2], R_sim, R_sim_V.T[-2])
    print >> out, 'expected rate:', -np.dot(v, np.diag(R))
    print >> out, 'second order expected rate:', -np.dot(v2, np.diag(R))
    print >> out
    print >> out, 'eigenvectors of R from eigenvectors of the similar matrix:'
    print >> out, R_sim_W
    print >> out, R_V_rebuilt
    print >> out
    return out.getvalue().rstrip()

def get_response_content(fs):
    return process(fs) + '\n'

