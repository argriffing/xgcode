"""For weighted Fiedler, compare node pseudo-duplication to explicit weights.

Build evidence for a weighted extension of the Fiedler cut theorem.
Compute the Fiedler valuation of nodes of a tree
augmented with a redundant tightly-connected node.
Compare this valuation to that of a valuation which uses weights explicitly.
As the strength of the pairing of the redundant nodes goes to infinity,
the valuations should become identical.
"""


from StringIO import StringIO
import math

import numpy as np

from SnippetUtil import HandlingError
import Form
import FormOut
import Euclid
import EigUtil
import MatrixUtil
import Util

g_D = np.array([
    [0, 1, 2, 3],
    [1, 0, 3, 4],
    [2, 3, 0, 5],
    [3, 4, 5, 0]], dtype=float)

def get_form():
    """
    @return: a list of form objects
    """
    form_objects = [
            Form.Matrix('D', 'distance matrix',
                g_D, MatrixUtil.assert_predistance),
            Form.Float('strength', 'pseudo-duplication connection strength',
                '1e6')]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_pseudoduplicate_laplacian(L, strength):
    """
    @param L: the laplacian matrix
    @param strength: the strength of the connection of the redundant pair
    @return: a new matrix with one more row and column
    """
    n = len(L)
    P = np.zeros((n+1, n+1))
    for i in range(n+1):
        for j in range(n+1):
            if i==n-1 and j==n-1:
                P[i,j] = L[i,j] + strength
            elif i==n and j==n:
                P[i,j] = strength
            elif (i==n and j==n-1) or (j==n and i==n-1):
                P[i,j] = -strength
            elif i < n and j < n:
                P[i,j] = L[i,j]
    return P

def get_response_content(fs):
    np.set_printoptions(linewidth=200)
    n = len(fs.D)
    # create the Laplacian matrix
    L = Euclid.edm_to_laplacian(fs.D)
    # create the Laplacian matrix with the extra node added
    L_dup = get_pseudoduplicate_laplacian(L, fs.strength)
    # get the principal axis projection from the Laplacian dup matrix
    X_w, X_v = EigUtil.principal_eigh(np.linalg.pinv(L_dup))
    L_dup_x = X_v * math.sqrt(X_w)
    # get masses summing to one
    m = np.array([1]*(n-1) + [2], dtype=float) / (n+1)
    # get the principal axis projection using the weight formula
    M = np.diag(np.sqrt(m))
    L_pinv = np.linalg.pinv(L)
    I = np.eye(n, dtype=float)
    E = I - np.outer(np.ones(n, dtype=float), m)
    ME = np.dot(M, E)
    Q = np.dot(ME, np.dot(L_pinv, ME.T))
    Q_w, Q_v = EigUtil.principal_eigh(Q)
    Q_x = Q_v * math.sqrt(Q_w) / np.sqrt(m)
    # make the response
    out = StringIO()
    print >> out, 'Laplacian matrix with pseudo-duplicate node:'
    print >> out, L_dup
    print >> out
    print >> out, 'principal axis projection:'
    print >> out, L_dup_x
    print >> out
    print >> out, 'principal axis projection using the weight formula:'
    print >> out, Q_x
    return out.getvalue()
