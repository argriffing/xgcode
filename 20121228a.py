"""
Check eigendecompositions of block triangular rate matrices.
"""

from StringIO import StringIO

import numpy as np
import scipy.linalg

import Form
import FormOut

def get_form():
    """
    @return: the body of a form
    """
    return [
            ]

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):
    np.set_printoptions(
            linewidth=1000000,
            threshold=1000000,
            )
    pre_Q = np.array([
        [0, 1, 1, 0],
        [1, 0, 0, 1],
        [0, 0, 0, 1],
        [0, 0, 1, 0],
        ], dtype=float)
    Q = pre_Q - np.diag(np.sum(pre_Q, axis=1))
    w, vl, vr = scipy.linalg.eig(Q, left=True, right=True)
    vl_inv = scipy.linalg.inv(vl)
    vr_inv = scipy.linalg.inv(vr)
    #
    out = StringIO()
    print >> out, 'Q:'
    print >> out, Q
    print >> out
    print >> out, 'w:'
    print >> out, w
    print >> out
    print >> out, 'vl:'
    print >> out, vl
    print >> out
    print >> out, 'vl.T:'
    print >> out, vl.T
    print >> out
    print >> out, 'inv(vl):'
    print >> out, vl_inv
    print >> out
    print >> out, 'vr:'
    print >> out, vr
    print >> out
    print >> out, 'vr.T:'
    print >> out, vr.T
    print >> out
    print >> out, 'inv(vr):'
    print >> out, vr_inv
    print >> out
    print >> out, 'inv(vl).T w vl.T:'
    print >> out, np.dot(vl_inv.T, np.dot(np.diag(w), vl.T))
    print >> out
    print >> out, 'vr w inv(vr):'
    print >> out, np.dot(vr, np.dot(np.diag(w), vr_inv))
    print >> out
    return out.getvalue()

