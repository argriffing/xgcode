"""
Check properties of a rate matrix augmented with hidden binary variables.

The binary variables allow or disallow transitions to corresponding states
of non-hidden variables.
For example, this could be used when some states are invisibly
disallowed at various times for whatever reason,
and the set of disallowed states changes
according to a continuous time process,
and only one state changes allowance status at any instant.
If the original rate matrix is Q with N states,
then the rate matrix augmented with these hidden binary variables
will be called Q* and it will have N * 2^(N-1) states.
The coefficient N corresponds to the N states of the visible variable,
and the term 2^(N-1) corresponds to the possible binary allowed/disallowed
states of the binary variables corresponding to the other N-1 states
of the visible variable.
The exponent is N-1 instead of N because the current state of the
visible variable cannot be disallowed.
So at most N-1 states of the visible variable can be disallowed
at any given time.
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
            Form.Float('blink_birth', 'tolerance birth rate', '0.01'),
            Form.Float('blink_death', 'tolerance death rate', '0.02'),
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
    if blink_birth < 0:
        raise Exception
    if blink_death < 0:
        raise Exception

    # define the original rate matrix
    #pre_Q = np.exp(np.random.randn(3, 3))
    #pre_Q = sample_reversible_pre_Q(4)
    """
    pre_Q = np.array([
        [0.0, 0.5],
        [0.1, 0.0],
        ], dtype=float)
    pre_Q = np.array([
        [0.0, 0.3, 0.1],
        [0.1, 0.0, 0.3],
        [0.1, 0.1, 0.0],
        ], dtype=float)
    pre_Q = np.array([
        [0.0, 3.0, 0.1],
        [0.1, 0.0, 3.0],
        [0.1, 0.1, 0.0],
        ], dtype=float)
    """
    pre_Q = np.array([
        [0.00, 6.00, 0.00],
        [0.00, 0.00, 6.00],
        [0.01, 0.00, 0.00],
        ], dtype=float)

    # construct the derived rate matrix
    pre_blink = binarytolerance.blinkize(pre_Q, blink_birth, blink_death)

    # check stationary distributions
    Q = binarytolerance.pre_Q_to_Q(pre_Q)
    W, V = scipy.linalg.eig(Q.T)
    v = V[:, np.argmin(np.abs(W))]
    v /= np.sum(v)
    print >> out, Q
    print >> out, W
    print >> out, v
    print >> out

    Q_blink = binarytolerance.pre_Q_to_Q(pre_blink)
    W, V = scipy.linalg.eig(Q_blink.T)
    v = V[:, np.argmin(np.abs(W))]
    v /= np.sum(v)
    print >> out, Q_blink
    print >> out, W
    print >> out, v
    print >> out

    n = pre_Q.shape[0]
    expo = binarytolerance.get_expanded_states(n)
    x = np.zeros(n)
    for i, (perm, mask) in enumerate(expo):
        x[perm[0]] += v[i]
    print >> out, x
    print >> out

    # show the result
    return out.getvalue()

