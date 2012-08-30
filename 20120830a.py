"""
Check an example of spectral tree inference by Zhang et al. 2011.
"""

from StringIO import StringIO
import math
from math import exp

import numpy as np
from scipy import linalg

import Form
import FormOut
import MatrixUtil
import StatsUtil
import Util

def get_form():
    """
    @return: the body of a form
    """
    return [
            ]

def get_form_out():
    return FormOut.Report()

def get_example_distance_matrix():
    """
    This is from a long branch attraction tree.
    In the true tree the first and last taxa are grouped together
    and are split from the middle two taxa.
    The idea is that bad algorithms will put the first two taxa together.
    """
    a = 9
    b = 5
    c = 2
    d = 1
    e = 1
    tri = np.array([
        [0,     a+e+b,  a+e+c,  a+d     ],
        [0,     0,      b+c,    b+e+d   ],
        [0,     0,      0,      c+e+d   ],
        [0,     0,      0,      0       ]], dtype=float)
    return tri + tri.T

def get_response_content(fs):
    D = get_example_distance_matrix()
    L_zhang = np.diag(np.sum(D, axis=1)) - D
    #
    out = StringIO()
    w, v = linalg.eigh(L_zhang)
    print >> out, D
    print >> out, w
    print >> out, v
    print >> out
    D_small = D[1:, 1:]
    L_small = np.diag(np.sum(D_small, axis=1)) - D_small
    w, v = linalg.eigh(L_small)
    print >> out, D_small
    print >> out, w
    print >> out, v
    print >> out
    return out.getvalue()

