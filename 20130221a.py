"""
Check more properties of a rate matrix augmented with hidden binary variables.

This assumes a non-reversible 3-state continuous time Markov process.
"""

from StringIO import StringIO
import itertools
import math

import numpy as np
import scipy.linalg

import Form
import FormOut
import MatrixUtil
from MatrixUtil import ndot
import combobreaker
import binarytolerance


def get_form():
    """
    @return: the body of a form
    """
    return [
            Form.Float('rate_01', 'rate A --> B',
                '1', low_exclusive=0),
            Form.Float('rate_12', 'rate B --> C',
                '0.01', low_exclusive=0),
            Form.Float('rate_20', 'rate C --> A',
                '0.01', low_exclusive=0),
            Form.Float('blink_birth', 'tolerance birth rate',
                '0.00000001', low_exclusive=0),
            Form.Float('blink_death', 'tolerance death rate',
                '0.0001', low_inclusive=0),
            ]

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):

    # define the amount of time we will search
    nseconds = 5

    # set up print options
    np.set_printoptions(
            linewidth=1000000,
            threshold=1000000,
            )
    out = StringIO()

    # set up some process parameters
    blink_birth = fs.blink_birth
    blink_death = fs.blink_death
    rate_01 = fs.rate_01
    rate_12 = fs.rate_12
    rate_20 = fs.rate_20

    # define the original rate matrix
    pre_Q = np.array([
        [0, rate_01, 0],
        [0, 0, rate_12],
        [rate_20, 0, 0],
        ], dtype=float)

    # construct the derived rate matrix
    pre_blink = binarytolerance.blinkize(pre_Q, blink_birth, blink_death)

    # check stationary distributions
    Q = binarytolerance.pre_Q_to_Q(pre_Q)
    W, V = scipy.linalg.eig(Q.T)
    v = V[:, np.argmin(np.abs(W))]
    v /= np.sum(v)
    #print >> out, Q
    #print >> out, W
    print >> out, 'original process stationary distribution:'
    print >> out, v
    print >> out

    Q_blink = binarytolerance.pre_Q_to_Q(pre_blink)
    W, V = scipy.linalg.eig(Q_blink.T)
    v = V[:, np.argmin(np.abs(W))]
    v /= np.sum(v)
    #print >> out, Q_blink
    #print >> out, W
    #print >> out, v
    #print >> out

    n = pre_Q.shape[0]
    expo = binarytolerance.get_expanded_states(n)
    x = np.zeros(n)
    for i, (perm, mask) in enumerate(expo):
        x[perm[0]] += v[i]
    print >> out, 'blink-ized stationary distribution:'
    print >> out, x
    print >> out

    # show the result
    return out.getvalue()

