"""
Check a matrix logarithm.
"""

from StringIO import StringIO
import itertools
import math

import numpy as np
import scipy.linalg

import Form
import FormOut
from MatrixUtil import ndot

def get_form():
    """
    @return: the body of a form
    """
    return [
            #Form.Integer('N', 'population size', 3, low=1, high=5),
            ]

def get_form_out():
    return FormOut.Report()

def hamdist(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)

def get_jeff_mut_trans(ascii_states, mu):
    k = len(ascii_states)
    mmu = np.zeros((k, k))
    for i, si in enumerate(ascii_states):
        for j, sj in enumerate(ascii_states):
            h = hamdist(si, sj)
            mmu[i, j] = (mu ** h) * ((1 - mu) ** (2 - h))
    return mmu

def get_mut_rate_matrix(ascii_states):
    k = len(ascii_states)
    pre_Q = np.zeros((k, k))
    for i, si in enumerate(ascii_states):
        for j, sj in enumerate(ascii_states):
            if hamdist(si, sj) == 1:
                pre_Q[i, j] = 1
    Q = pre_Q - np.diag(np.sum(pre_Q, axis=1))
    return Q

def get_equivalent_rate(mu):
    """
    Get a rate that gives an expm matrix with the same entries as Jeff has.
    """
    return -0.25 * math.log(1 - 4*mu*(1-mu))

def get_response_content(fs):

    mu = 0.025
    ascii_states = ['AB', 'Ab', 'aB', 'ab']
    k = len(ascii_states)

    # set up print options
    np.set_printoptions(
            linewidth=1000000,
            threshold=1000000,
            )
    out = StringIO()

    # show the original transition matrix
    mmu = get_jeff_mut_trans(ascii_states, mu)
    print >> out, 'original transition matrix:'
    print >> out, mmu
    print >> out
    print >> out, 'logm of original transition matrix:'
    print >> out, scipy.linalg.logm(mmu)
    print >> out
    print >> out

    # show a transition matrix derived from a rate matrix
    Q = get_mut_rate_matrix(ascii_states)
    print >> out, 'a scaled rate matrix:'
    print >> out, mu*Q
    print >> out
    print >> out, 'expm of the scaled rate matrix:'
    print >> out, scipy.linalg.expm(mu*Q)
    print >> out
    print >> out

    # do a thing
    mu_adjusted = get_equivalent_rate(mu)
    print >> out, 'adjusted mu:'
    print >> out, mu_adjusted
    print >> out
    print >> out, 'expm of the scaled adjusted rate matrix:'
    print >> out, scipy.linalg.expm(mu_adjusted*Q)
    print >> out

    # show the result
    return out.getvalue()

