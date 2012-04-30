"""
Look at relaxation times of aggregated independent site reversible CTMC models.

Two independently evolving sites are aggregated.
"""

from StringIO import StringIO
import math
import itertools
import random

import numpy as np
import scipy
from scipy import linalg

import Form
import FormOut
from MatrixUtil import ndot
import mrate
import divtime

def get_form():
    """
    @return: the body of a form
    """
    # define the form objects
    form_objects = [
            Form.Integer('nstates', 'use this many states per site',
                3, low=2, high=4)]
    return form_objects

def get_form_out():
    return FormOut.Report()

def aggregate(R):
    n = len(R)
    Q = np.zeros((n*n, n*n))
    for a in range(n):
        for b in range(n):
            for na in range(n):
                if na != a:
                    Q[a * n + b, na * n + b] = R[a, na]
            for nb in range(n):
                if nb != b:
                    Q[a * n + b, a * n + nb] = R[b, nb]
    Q -= np.diag(np.sum(Q, axis=1))
    return Q

def get_response_content(fs):
    out = StringIO()
    np.set_printoptions(linewidth=200)
    # get the user defined variables
    n = fs.nstates
    # sample a random reversible CTMC rate matrix
    v = divtime.sample_distribution(n)
    S = divtime.sample_symmetric_rate_matrix(n)
    R = mrate.to_gtr_halpern_bruno(S, v)
    distn = mrate.R_to_distn(R)
    spectrum = scipy.linalg.eigvalsh(mrate.symmetrized(R))
    print >> out, 'random reversible CTMC rate matrix:'
    print >> out, R
    print >> out
    print >> out, 'stationary distribution:'
    print >> out, distn
    print >> out
    print >> out, 'spectrum:'
    print >> out, spectrum
    print >> out
    Q = aggregate(R)
    distn = mrate.R_to_distn(Q)
    spectrum = scipy.linalg.eigvalsh(mrate.symmetrized(Q))
    print >> out, 'aggregated rate matrix:'
    print >> out, Q
    print >> out
    print >> out, 'stationary distribution:'
    print >> out, distn
    print >> out
    print >> out, 'spectrum:'
    print >> out, spectrum
    print >> out
    return out.getvalue()

