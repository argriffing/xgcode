"""
This module is about general finite-state continuous-time Markov processes.

The python packages numpy and scipy are used throughout.
At some point I should probably make separate modules
to emphasize the distinction between reversible and general
continuous-time Markov processes.
"""

import math
import unittest

import numpy as np
import scipy
from scipy import linalg

from MatrixUtil import ndot

def expm_spectral(R, t):
    """
    This is for testing expm_diff_spectral only.
    You should use scipy.linalg.expm instead.
    """
    n = len(R)
    v = R_to_distn(R)
    S = symmetrized(R)
    w, U = np.linalg.eigh(S)
    P = np.zeros_like(R)
    for i in range(n):
        for j in range(n):
            for k in range(n):
                a = (v[j] / v[i])**0.5
                b = U[i, k] * U[j, k]
                c = math.exp(t * w[k])
                P[i, j] += a * b * c
    return P

def expm_diff_spectral(R, t):
    """
    Get the rates of change of transition probabilities at time t.
    Use the spectral representation.
    @return: entrywise derivative of transition matrix at time t
    """
    n = len(R)
    v = R_to_distn(R)
    S = symmetrized(R)
    w, U = np.linalg.eigh(S)
    P_diff = np.zeros_like(R)
    for i in range(n):
        for j in range(n):
            for k in range(n):
                a = (v[j] / v[i])**0.5
                b = U[i, k] * U[j, k]
                c = w[k] * math.exp(t * w[k])
                P_diff[i, j] += a * b * c
    return P_diff

def _R_to_eigenpair(R):
    n = len(R)
    Wl, Vl = scipy.linalg.eig(R, left=True, right=False)
    val_vec_pairs = [(abs(Wl[i]), Vl[:,i]) for i in range(n)]
    r_recip, fiedler = sorted(val_vec_pairs)[1]
    return r_recip, fiedler

def R_to_fiedler(R):
    r_recip, fiedler = _R_to_eigenpair(R)
    return fiedler

def R_to_relaxation_time(R):
    r_recip, fiedler = _R_to_eigenpair(R)
    return 1 / r_recip

def R_to_distn(R):
    """
    @param R: rate matrix
    @return: stationary distribution
    """
    n = len(R)
    Wl, Vl = scipy.linalg.eig(R, left=True, right=False)
    val_vec_pairs = [(abs(Wl[i]), Vl[:,i]) for i in range(n)]
    dummy, pi_eigenvector = min(val_vec_pairs)
    total = np.sum(pi_eigenvector)
    pi_arr = np.array([v/total for v in pi_eigenvector])
    return pi_arr

def R_to_total_rate(R):
    #TODO: the jargon should be 'expected rate' rather than 'total rate'
    n = len(R)
    distn = R_to_distn(R)
    total_rate = 0.0
    for i in range(n):
        total_rate -= distn[i] * R[i, i]
    return total_rate

def symmetrized(R):
    """
    Get the symmetrized matrix of a reversible markov process.
    This returns a symmetric matrix that is not a rate matrix
    because rows do not sum to zero.
    """
    v = R_to_distn(R)
    lam = np.diag(np.sqrt(v))
    rlam = np.diag(np.reciprocal(np.sqrt(v)))
    return ndot(lam, R, rlam)

class TestMrate(unittest.TestCase):

    def test_expm(self):
        M = np.array([
            [-1.55273124,  0.40323905,  0.90129456,  0.24819763],
            [ 2.41191856, -5.14753325,  2.29550348,  0.44011122],
            [ 0.82711917,  0.3521918,  -1.44057109,  0.26126012],
            [ 3.01693453,  0.89439764,  3.46050956, -7.37184173]])
        t = 0.3
        observed = scipy.linalg.expm(M * t)
        expected = expm_spectral(M, t)
        self.assertTrue(np.allclose(observed, expected))

if __name__ == '__main__':
    unittest.main()
