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

def get_scaled_eigenvectors(C_full, diploid_and_biallelic):
    """
    Scale using the eigenvalues.
    @param C_full: matrix of float counts where each row represents an OTU
    @param diploid_and_biallelic: a flag
    @return: eigenvectors scaled by eigenvalues
    """
    evals, evecs = get_eval_evec_pairs(C_full, diploid_and_biallelic)
    return [w*v for w, v in zip(evals, evecs)]

def get_sqrt_scaled_eigenvectors(C_full, diploid_and_biallelic):
    """
    Scale using the square roots of eigenvalues.
    @param C_full: matrix of float counts where each row represents an OTU
    @param diploid_and_biallelic: a flag
    @return: eigenvectors scaled by eigenvalues
    """
    evals, evecs = get_eval_evec_pairs(C_full, diploid_and_biallelic)
    # Eigenvalues are from a gramian matrix
    # so they should theoretically be all nonnegative.
    # From the numerical computation, some might be negative and near zero.
    # Therefore change negative eigenvalues to zero.
    evals[np.where(evals<0)] = 0.0
    return [math.sqrt(w)*v for w, v in zip(evals, evecs)]

def get_unscaled_eigenvectors(C_full, diploid_and_biallelic):
    """
    Get the unscaled eigenvectors.
    That is, the eigenvectors will be orthonormal.
    @param C_full: matrix of float counts where each row represents an OTU
    @param diploid_and_biallelic: a flag
    @return: unscaled eigenvectors
    """
    return evecs

def get_eval_evec_pairs(C_full, diploid_and_biallelic):
    """
    Input rows are OTUs and columns are loci.
    Each element of the input data is a count.
    @param C_full: matrix of float counts where each row represents an OTU
    @param diploid_and_biallelic: a flag
    @return: (eigenvalues, eigenvectors)
    """
    # create the floating point count matrix
    m_full, n_full = C_full.shape
    # check compatibility of counts and ploidy
    if diploid_and_biallelic:
        if np.max(C_full) > 2:
            raise ValueError(
                    'no count should be greater than two for diploid data')
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
        variances = p * (1 - p)
        M /= np.sqrt(variances)
    # construct the sample covariance matrix
    # FIXME this should probably use a singular value decomposition instead
    X = np.dot(M, M.T) / n
    # get the eigendecomposition of the covariance matrix
    return EigUtil.eigh(X)
