"""
Approximate spectral gap change with respect to the stationary distribution.

Changes in stationary distribution are assumed to affect
the rate matrix and hence the spectral gap
through the f=1/2 approximation of the Sella-Hirsh approximation
of natural selection.
This turned out to work better than I had expected.
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
            Form.Float('eps', 'norm of stationary distribution change',
                0.001, low_exclusive=0, high_exclusive=1),
            Form.RadioGroup('seltype', 'selection type', [
                Form.RadioItem('knudsen', 'Knudsen-Miyamoto', True),
                Form.RadioItem('sella', 'Sella-Hirsh')]),
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
    # Sample some numbers then subtract mean then normalize.
    dv = np.random.exponential(1, n)
    dv -= np.mean(dv)
    dv *= fs.eps / np.dot(dv, dv)
    qv = v + dv
    if any(qv < 0) or any(1 < qv):
        raise ValueError(
            'the stationary distribution change was too large '
            'for the randomly sampled process')
    qpsi = np.sqrt(qv)
    # define the rate matrix
    if fs.knudsen:
        Q = (S.T / qpsi).T * qpsi
    elif fs.sella:
        Q = R.copy()
        for a in range(n):
            for b in range(n):
                if a != b:
                    tau = (qv[b] / v[b]) / (qv[a] / v[a])
                    Q[a, b] *= math.log(tau) / (1 - 1/tau)
    Q -= np.diag(np.sum(Q, axis=1))
    # construct the symmetric matrix that is similar to Q
    Q_sim = (Q.T * qpsi).T / qpsi
    Q_sim_W, Q_sim_V = scipy.linalg.eigh(Q_sim)
    Q_gap = -Q_sim_W[-2]
    # report some stuff
    print >> out, 'sampled rate matrix R:'
    print >> out, R
    print >> out, 'stationary distn:', v
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
    print >> out, 'mutation-selection balance matrix Q:'
    print >> out, Q
    print >> out, 'stationary distn:', qv
    print >> out, 'spectral gap:', Q_gap
    print >> out
    print >> out, 'symmetric matrix similar to Q:'
    print >> out, Q_sim
    print >> out
    print >> out, 'pi(Q) - pi(R):', dv
    print >> out, 'gap(Q) - gap(R):', Q_gap - R_gap
    print >> out, 'diag(Q) - diag(R):', np.diag(Q) - np.diag(R)
    print >> out, 'trace(Q) - trace(R):', np.trace(Q) - np.trace(R)
    print >> out
    print >> out, 'rate away estimate of spectral gap change:'
    print >> out, np.dot(np.diag(Q) - np.diag(R), R_sim_V.T[-2]**2)
    print >> out
    return out.getvalue().rstrip()

def get_response_content(fs):
    return process(fs) + '\n'

