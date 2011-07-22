"""
Look at a few matrices associated with a path graph.
"""

from StringIO import StringIO
import numpy as np
import scipy
import scipy.linalg

import Form
import FormOut
import MatrixUtil

def get_form():
    """
    @return: the body of a form
    """
    # define the form objects
    form_objects = [
            Form.Float('edge_length', 'length per edge',
                2.0, low_exclusive=0),
            Form.Integer('nvertices', 'number of vertices',
                4, low=2, high=10)]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):
    # make the laplacian matrix for the graph
    weight = 1 / fs.edge_length
    n = fs.nvertices
    L = np.zeros((n,n))
    # set the diagonal
    for i in range(n):
        L[i,i] = 2 * weight
    L[0,0] = weight
    L[-1,-1] = weight
    # complete the tridiagonal
    for i in range(n-1):
        L[i+1,i] = -weight
        L[i,i+1] = -weight
    # define other matrices
    L_pinv = np.linalg.pinv(L)
    HDH = -2*L_pinv
    v = np.diag(HDH)
    e = np.ones(n)
    D = HDH - (np.outer(v, e) + np.outer(e, v))/2
    # show some matrices
    out = StringIO()
    np.set_printoptions(linewidth=300)
    print >> out, 'Laplacian matrix:'
    print >> out, L
    print >> out
    print >> out, 'HDH:'
    print >> out, HDH
    print >> out
    print >> out, 'EDM:'
    print >> out, D
    print >> out
    return out.getvalue()

