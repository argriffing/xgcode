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

def two_state_expm(a, b, t):
    """
    Compute the matrix exponential of a continuous time Markov process.
    The Markov process has only two states and is defined by the two rates.
    @param a: rate of change from first state to second state
    @param b: rate of change from second state to first state
    @param t: elapsed time
    @return: a conditional transition matrix
    """
    r = a + b
    p = scipy.exp(-r*t)
    P = np.array([
        [p*a + b, a - a*p],
        [b - b*p, a + b*p],
        ], dtype=float) / r
    return P

def get_pq_init_off(a, w, r, t):
    """
    Compute differential equations evaluated at the given time.
    @param a: rate from off to on
    @param w: rate from on to off
    @param r: poisson event rate
    @param t: elapsed time
    @return: p, q
    """
    x = scipy.sqrt(scipy.square(a + r + w) - 4*a*r)
    denom = 2 * x * scipy.exp(t * (x + a + r + w) / 2)
    p = (scipy.exp(t*x)*(x + r + w - a) + (x - r - w + a)) / denom
    q = (2 * a * (scipy.exp(t * x) - 1)) / denom
    return p, q

def get_pq_init_on(a, w, r, t):
    """
    Compute differential equations evaluated at the given time.
    @param a: rate from off to on
    @param w: rate from on to off
    @param r: poisson event rate
    @param t: elapsed time
    @return: p, q
    """
    x = scipy.sqrt(scipy.square(a + r + w) - 4*a*r)
    denom = 2 * x * scipy.exp(t * (x + a + r + w) / 2)
    p = (2 * w * (scipy.exp(t * x) - 1)) / denom
    q = (scipy.exp(t*x)*(x - r - w + a) + (x + r + w - a)) / denom
    return p, q

def expm_shortcut(a, w, r, t):
    """
    This is a simpler approach that uses matrix exponential.
    It uses a continuous time Markov process with three states.
    @param a: rate from off to on
    @param w: rate from on to off
    @param r: poisson event rate
    @param t: elapsed time
    """
    pre_Q = np.array([
        [0, a, 0],
        [w, 0, r],
        [0, 0, 0],
        ], dtype=float)
    Q = pre_Q - np.diag(np.sum(pre_Q, axis=1))
    P = scipy.linalg.expm(Q*t)
    return P

def get_response_content(fs):

    # set up print options
    np.set_printoptions(
            linewidth=1000000,
            threshold=1000000,
            )
    out = StringIO()

    # read the specified parameter values
    a = fs.lam
    w = fs.mu
    r = fs.alpha
    t = fs.duration
    if fs.initial_on:
        b_0 = 1
    if fs.initial_off:
        b_0 = 0
    if fs.final_on:
        b_t = 1
    if fs.final_off:
        b_t = 0

    """
    # You should be able to prove that this value
    # is real when all of a, r, w are positive.
    x_sqr = scipy.square(a + r + w) - 4*a*r
    x = scipy.sqrt(x_sqr)

    # Write the equivalent expression for the square of x
    # in a way that is more obviously positive.
    x_sqr_b = (scipy.square(scipy.sqrt(a) - scipy.sqrt(r)) + w) * (
            a + r + w + 2*scipy.sqrt(a*r))

    # Compute the solutions of differential equations, using Mathematica.
    # The corresponding system of differential equations is
    # p'[t] == -a*p[t] + w*q[t]
    # q'[t] == a*p[t] - (w+r)*q[t]
    # p[0] == 1
    # q[0] == 0
    #
    denom = 2 * x * scipy.exp(t * (x + a + r + w) / 2)
    p = (scipy.exp(t*x)*(x + r + w - a) + (x - r - w + a)) / denom
    q = (2 * a * (scipy.exp(t * x) - 1)) / denom
    """

    if b_0 == 1:
        p, q = get_pq_init_on(a, w, r, t)
    if b_0 == 0:
        p, q = get_pq_init_off(a, w, r, t)

    # remember that
    # p is no trans and off
    # q is no trans and on
    conditional_markov = two_state_expm(a, w, t)[b_0, b_t]
    if b_t == 0:
        endpoint_conditioned_prob = p / conditional_markov
    if b_t == 1:
        endpoint_conditioned_prob = q / conditional_markov
    
    # report some values
    print >> out, 'differential equation solutions at the given time'
    print >> out, 'p(t):', p
    print >> out, 'q(t):', q
    print >> out
    #print >> out, 'to prove no imaginary numbers, these two things are equal'
    #print >> out, x_sqr
    #print >> out, x_sqr_b
    #print >> out
    print >> out, 'probability of no Poisson event unconditional on final state'
    print >> out, 'p(t) + q(t):', p + q
    print >> out
    print >> out, 'endpoint conditioned prob of no Poisson event'
    print >> out, endpoint_conditioned_prob
    print >> out
    print >> out, 'matrix exponential shortcut:'
    print >> out, expm_shortcut(a, w, r, t)
    print >> out

    # show the result
    return out.getvalue()

