"""
Friendlier eigendecomposition.
This module is unnecessary because the linalg module of the scipy package
already does more intelligent things with eigendecomposition than numpy does.
"""

import numpy as np
import scipy
from scipy import linalg

def eigh(M):
    """
    Unlike numpy.eigh, the results of this function are ordered.
    The returned eigenvalues and eigenvectors are sorted
    by decreasing eigenvalue.
    The eigenvectors are returned as a list of one dimensional numpy arrays.
    If treated as a 2D array, the returned eigenvector list is the transpose
    of the eigenvector result of numpy.eig.
    The eigenvectors are normalized unit right eigenvectors.
    @param M: a symmetric numpy 2D array
    @return: eigenvalues and eigenvectors
    """
    # FIXME the docstring describes the transpose of the returned vectors
    """
    w, vt = np.linalg.eigh(M)
    ws, vs = w.tolist(), vt.T.tolist()
    sorted_pairs = list(reversed(sorted(zip(ws, vs))))
    w, v = zip(*sorted_pairs)
    return np.array(w), [np.array(x) for x in v]
    """
    #return scipy.linalg.eigh(M)
    W, VT = scipy.linalg.eigh(M)
    W_reversed = np.array(list(reversed(W)))
    V_reversed = np.array(list(reversed(VT.T)))
    return W_reversed, V_reversed

def principal_eigh(M):
    """
    @param M: a symmetric numpy 2D array
    @return: principal eigenvalue and eigenvector
    """
    #
    """
    w, vt = np.linalg.eigh(M)
    ws, vs = w.tolist(), vt.T.tolist()
    val, vec = max(zip(ws, vs))
    return val, np.array(vec)
    """
    W, VT = scipy.linalg.eigh(M)
    return W[-1], VT.T[-1]

