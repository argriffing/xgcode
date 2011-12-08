import math

import numpy as np
import scipy
from scipy import linalg

from MatrixUtil import ndot

def R_to_relaxation_time(R):
    n = len(R)
    Wl, Vl = scipy.linalg.eig(R, left=True, right=False)
    val_vec_pairs = [(abs(Wl[i]), Vl[:,i]) for i in range(n)]
    r_recip, fiedler = sorted(val_vec_pairs)[1]
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
    n = len(R)
    distn = R_to_distn(R)
    total_rate = 0.0
    for i in range(n):
        total_rate -= distn[i] * R[i, i]
    return total_rate

def symmetrized(R):
    """
    Get the symmetrized matrix.
    This returns a symmetric matrix that is not a rate matrix
    because rows do not sum to zero.
    """
    v = R_to_distn(R)
    lam = np.diag(np.sqrt(v))
    rlam = np.diag(np.reciprocal(np.sqrt(v)))
    return ndot(lam, R, rlam)

