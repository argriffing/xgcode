"""For weighted Fiedler, compare exact node duplication to explicit weights.

Build evidence for a weighted extension of the Fiedler cut theorem.
Compute the Fiedler valuation of nodes of a tree
augmented with a node with distance zero from an existing node.
Compare this valuation to that of a valuation which uses weights explicitly.
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
                g_D, MatrixUtil.assert_predistance)]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_duplicate_edm(D):
    """
    @param D: the distance matrix
    @return: a new distance matrix with one more row and column
    """
    n = len(D)
    P = np.zeros((n+1, n+1))
    for i in range(n+1):
        for j in range(n+1):
            if i==n and j==n:
                continue
            elif i==n:
                P[i,j] = D[i-1,j]
            elif j==n:
                P[i,j] = D[i,j-1]
            else:
                P[i,j] = D[i,j]
    return P

def get_response_content(fs):
    np.set_printoptions(linewidth=200)
    n = len(fs.D)
    # create the distance matrix with the extra node added
    D_dup = get_duplicate_edm(fs.D)
    # get the principal axis projection from the distance dup matrix
    HSH = -(0.5) * MatrixUtil.double_centered(D_dup)
    X_w, X_v = EigUtil.principal_eigh(HSH)
    D_dup_x = X_v * math.sqrt(X_w)
    # get masses summing to one
    m = np.array([1]*(n-1) + [2], dtype=float) / (n+1)
    # get the principal axis projection using the weight formula
    M = np.diag(np.sqrt(m))
    I = np.eye(n, dtype=float)
    E = I - np.outer(np.ones(n, dtype=float), m)
    ME = np.dot(M, E)
    Q = -(0.5) * np.dot(ME, np.dot(fs.D, ME.T))
    Q_w, Q_v = EigUtil.principal_eigh(Q)
    Q_x = Q_v * math.sqrt(Q_w) / np.sqrt(m)
    # make the response
    out = StringIO()
    print >> out, 'distance matrix with exact duplicate node:'
    print >> out, D_dup
    print >> out
    print >> out, 'principal axis projection:'
    print >> out, D_dup_x
    print >> out
    print >> out, 'principal axis projection using the weight formula:'
    print >> out, Q_x
    return out.getvalue()
