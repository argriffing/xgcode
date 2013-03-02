"""
Check probability of zero Poisson events from a 2-state MMPP.

The 2-state Markov modulated poisson process (MMPP)
is defined by two Markov rate parameters and two Poisson event rates.
Note that a 2-state Markov process is time-reversible.
For our purposes, we will treat one of the two Markov states as active
and the other as inactive, and the inactive Markov state will have
zero rate of emitting a Poisson event.
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

# FullSimplify[
#   DSolve[ {
#     p'[t] == -a*p[t] + w*q[t],
#     q'[t] == a*p[t] - (w+r)*q[t],
#     p[0]==1, q[0]==0
#   }, {p[t], q[t]}, t]
# ]
#
#

def get_form():
    """
    @return: the body of a form
    """
    return [
            Form.Float('lam', 'rate: off -> on',
                0.2, low_exclusive=0),
            Form.Float('mu', 'rate: on -> off',
                0.1, low_exclusive=0),
            Form.Float('alpha', 'active state Poisson rate',
                0.5, low_exclusive=0),
            Form.Float('duration', 'duration',
                2.0, low_exclusive=0),
            Form.RadioGroup('initial_activation', 'initial activation state', [
                Form.RadioItem('initial_on', 'on', True),
                Form.RadioItem('initial_off', 'off'),
                ]),
            Form.RadioGroup('final_activation', 'final activation state', [
                Form.RadioItem('final_on', 'on', True),
                Form.RadioItem('final_off', 'off'),
                ]),
            ]

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):

    # set up print options
    np.set_printoptions(
            linewidth=1000000,
            threshold=1000000,
            )
    out = StringIO()

    # read the specified parameter values
    a = fs.lam
    r = fs.alpha
    w = fs.mu
    t = fs.duration

    # You should be able to prove that this value
    # is real when all of a, r, w are positive.
    x_sqr = scipy.square(a + r + w) - 4*a*r
    x = scipy.sqrt(x_sqr)

    # Write the equivalent expression for x.
    x_sqr_b = (scipy.square(scipy.sqrt(a) - scipy.sqrt(r)) + w) * (
            a + r + w + 2*scipy.sqrt(a*r))

    # report some values
    denom = 2 * x * scipy.exp(t * (x + a + r + w) / 2)
    p = (scipy.exp(t*x)*(x + r + w - a) + (x - r - w + a)) / denom
    q = (2 * a * (scipy.exp(t * x) - 1)) / denom
    print >> out, 'differential equation solutions at the given time'
    print >> out, 'p(t):', p
    print >> out, 'q(t):', q
    print >> out
    print >> out, 'to prove no imaginary numbers, these two things are equal'
    print >> out, x_sqr
    print >> out, x_sqr_b
    print >> out

    # show the result
    return out.getvalue()

