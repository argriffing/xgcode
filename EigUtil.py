"""
Friendlier eigendecomposition.
"""

import numpy as np

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
    w, vt = np.linalg.eigh(M)
    ws, vs = w.tolist(), vt.T.tolist()
    sorted_pairs = list(reversed(sorted(zip(ws, vs))))
    w, v = zip(*sorted_pairs)
    return np.array(w), [np.array(x) for x in v]