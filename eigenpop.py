"""
Analyses of eigendecompositions of population genetic matrices.

This code is from the paper
Population Structure and Eigenanalysis
by Nick Patterson, Alkes Price, and David Reich.
In this code the word matrix really means 2D numpy array.
"""

import math

import numpy as np

import EigUtil

def get_scaled_eigenvectors(count_matrix, diploid_and_biallelic):
    """
    Input rows are OTUs and columns are loci.
    Each element of the input data is a count.
    @param count_matrix: matrix of counts where each row represents an OTU
    @param diploid_and_biallelic: a flag
    @return: eigenvectors scaled by square roots of eigenvalues
    """
    # create the floating point count matrix
    C_full = np.array(count_matrix, dtype=float)
    m_full, n_full = C_full.shape
    # check compatibility of counts and ploidy
    if diploid_and_biallelic:
        if np.max(C_full) > 2:
            msg = 'no count should be greater than two for diploid data'
            raise ValueError(msg)
    # remove invariant columns
    C = np.vstack([v for v in C_full.T if len(set(v))>1]).T
    # get the shape of the matrix
    m, n = C.shape
    # get the column means
    u = C.mean(axis=0)
    # get the centered and normalized counts matrix
    M = (C - u)
    # normalize if diploid and biallelic
    if diploid_and_biallelic:
        p = u/2
        M /= np.sqrt(p * (1 - p))
    # construct the sample covariance matrix
    X = np.dot(M, M.T) / n
    # get the eigendecomposition of the covariance matrix
    evals, evecs = EigUtil.eigh(X)
    # scale the eigenvectors by the square roots of the eigenvalues
    pcs = [math.sqrt(w)*v for w, v in zip(evals, evecs)]
    # return the eigenvectors
    return pcs
