"""Check a correlation matrix property.

Make a data matrix X and try to find a data matrix Y such
that the correlation matrix corresponding to Y is equal to the
elementwise squared correlation matrix corresponding to X.
"""

from StringIO import StringIO
import random
import math

import numpy as np

from SnippetUtil import HandlingError
import MatrixUtil
import Form
import FormOut
import NewickIO

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            #Form.CheckGroup('options', 'options', [
                #Form.CheckItem('standardize', 'standardize the variables to 1 variance', True)]),
            Form.Integer('n', 'n', 4, low=1, high=10),
            Form.Integer('p', 'p', 20, low=1, high=1000)]
    return form_objects

"""
    print >> out, 'U:'
    print >> out, U.shape
    print >> out, U
    print >> out
    print >> out, 'S:'
    print >> out, S.shape
    print >> out, S
    print >> out
    print >> out, 'VT:'
    print >> out, VT.shape
    print >> out, VT
    print >> out
"""

def get_form_out():
    return FormOut.Report()

def get_bilinear_form(M, v):
    """
    @param M: a correlation matrix
    @param v: a conformant vector
    """
    n = len(v)
    accum = 0
    for i in range(n):
        for j in range(n):
            accum += v[i] * v[j] * M[i, j]
    return accum

def get_dominant_ev_pair(M):
    """
    @param M: a correlation matrix
    @return: an eigenvalue eigenvector pair
    """
    w, vt = np.linalg.eigh(M)
    return max(zip(w, vt.T))

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    out = StringIO()
    # check the input
    n = fs.n
    p = fs.p
    if p < n * n:
        raise HandlingError('p should be at least as big as n^2')
    # make a data matrix with random values
    X_raw = np.zeros((p, n))
    for i in range(p):
        for j in range(n):
            X_raw[i,j] = random.expovariate(1.0)
    X = X_raw.copy()
    # force each row to have mean of zero
    for i in range(p):
        row_mean = np.mean(X[i])
        X[i] -= row_mean
    # force each row to have variance of one
    for i in range(p):
        row_std = np.std(X[i])
        X[i] /= row_std
    # account for the number of columns
    X /= math.sqrt(n)
    # create the covariance matrix associated with X (correlation if standardized)
    X_corr = np.dot(X, X.T)
    # assert that XX' is a the right matrix
    assert np.allclose(X_corr, np.corrcoef(X_raw))
    # create the Y matrix from the standardized X matrix
    Y = np.zeros((p, n*n))
    for i in range(p):
        for j in range(n):
            for k in range(n):
                Y[i,j*n + k] = X[i,j] * X[i,k]
    # create the correlation matrix associated with Y
    Y_corr = np.dot(Y, Y.T)
    # Assert that the correlation matrix associated with Y is
    # equal to the entrywise square of the correlation matrix associated with X.
    assert np.allclose(Y_corr, X_corr ** 2.0)
    # create the Z matrix from the standardized X matrix
    ncombos = (n*(n+1))/2
    Z = np.zeros((p, ncombos))
    for i in range(p):
        index = 0
        for j in range(n):
            Z[i,index] = X[i,j] ** 2
            index += 1
            for k in range(j):
                Z[i,index] = math.sqrt(2) * X[i,j] * X[i,k]
                index += 1
    # create the correlation matrix associated with Z
    Z_corr = np.dot(Z, Z.T)
    # Assert that the correlation matrix associated with Z is
    # equal to the entrywise square of the correlation matrix associated with X.
    assert np.allclose(Z_corr, X_corr ** 2.0)
    # create the W matrix from the standardized X matrix
    U, S_array, VT = np.linalg.svd(X, full_matrices=0)
    S = np.diag(S_array)
    assert np.allclose(X, np.dot(U, np.dot(S, VT)))
    X_reduced = np.array([row[:-1] for row in np.dot(U, S)])
    ncombos_reduced = (n*(n-1))/2
    W = np.zeros((p, ncombos_reduced))
    for i in range(p):
        index = 0
        for j in range(n-1):
            W[i,index] = X_reduced[i,j] ** 2
            index += 1
            for k in range(j):
                W[i,index] = math.sqrt(2) * X_reduced[i,j] * X_reduced[i,k]
                index += 1
    # create the correlation matrix associated with Z
    W_corr = np.dot(W, W.T)
    # Assert that the correlation matrix associated with Z is
    # equal to the entrywise square of the correlation matrix associated with X.
    assert np.allclose(W_corr, X_corr ** 2.0)
    # show some results
    print >> out, 'X:'
    print >> out, X
    print >> out
    print >> out, 'the correlation matrix for X:'
    print >> out, X_corr
    print >> out
    print >> out, 'entrywise square of the correlation matrix for X:'
    print >> out, X_corr ** 2.0
    print >> out
    print >> out, 'Y:'
    print >> out, Y
    print >> out
    print >> out, 'the correlation matrix for Y:'
    print >> out, np.dot(Y, Y.T)
    print >> out
    print >> out, 'Z:'
    print >> out, Z
    print >> out
    print >> out, 'W:'
    print >> out, W
    print >> out
    print >> out, 'eigenvalues of the correlation matrix for X:'
    print >> out, np.linalg.eigvalsh(X_corr)
    print >> out
    print >> out, 'eigenvalues of the correlation matrix for Y:'
    print >> out, np.linalg.eigvalsh(Y_corr)
    print >> out
    print >> out, 'dominant ev pair from X_corr:'
    print >> out, get_dominant_ev_pair(X_corr)
    print >> out
    print >> out, 'dominant ev pair from Y_corr:'
    print >> out, get_dominant_ev_pair(Y_corr)
    print >> out
    print >> out, 'dominant eigenvector from Y_corr applied to X_corr:'
    print >> out, get_bilinear_form(X_corr, get_dominant_ev_pair(Y_corr)[1])
    print >> out
    print >> out, 'dominant eigenvector from Y_corr applied to Y_corr:'
    print >> out, get_bilinear_form(Y_corr, get_dominant_ev_pair(Y_corr)[1])
    print >> out
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
