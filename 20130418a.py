"""
Check a diagonalization of a structured three-state rate matrix.
"""

from StringIO import StringIO
import random
import math

import numpy as np

import Form
import FormOut

def random_category(distn):
    """
    Sample from a categorical distribution.
    Note that this is not the same as random.choice(distn).
    Maybe a function like this will eventually appear
    in python or numpy or scipy.
    @param distn: categorical distribution as a stochastic vector
    @return: category index as a python integer
    """
    nstates = len(distn)
    np_index = np.dot(np.arange(nstates), np.random.multinomial(1, distn))
    return int(np_index)


def get_form():
    """
    @return: the body of a form
    """
    return [
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

    # sample the random rates
    a = random.expovariate(1)
    b = random.expovariate(1)
    c = random.expovariate(1)

    # rate matrix
    Q = np.array([
        [-a, a, 0],
        [b, -(b+c), c],
        [0, 0, 0],
        ], dtype=float)

    # Compute the square root of a thing that is non-negative
    # when all rates are non-negative.
    # The non-negativity could be proved by factoring
    # the difference of squares and then using the
    # inequality of arithmetic and geometric means.
    d = math.sqrt((a+b+c)*(a+b+c) - 4*a*c)

    # According to wolfram alpha this is the diagonalization Q = U diag(w) V
    # where U V = I.
    U = np.array([
        [1, (-a+b+c-d)/(2*b), (-a+b+c+d)/(2*b)],
        [1, 1, 1],
        [1, 0, 0],
        ], dtype=float)
    w = np.array([
        0,
        -(a+b+c+d)/2,
        -(a+b+c-d)/2,
        ], dtype=float)
    V = np.array([
        [0, 0, 1],
        [-b/d, (-a+b+c+d)/(2*d), (a+b-c-d)/(2*d)],
        [b/d, (a-b-c+d)/(2*d), (-a-b+c-d)/(2*d)],
        ], dtype=float)

    # now lets get some path samples

    # define the initial distribution
    initial_distn = np.array([
        Q[1, 0] / (Q[0, 1] + Q[1, 0]),
        Q[0, 1] / (Q[0, 1] + Q[1, 0]),
        ], dtype=float)
    
    # define distributions conditional on transitioning
    P = Q - np.diag(np.diag(Q))
    for i in range(P.shape[0]):
        P[i] /= np.sum(P[i])

    # initial state, final state, source state, sink state
    transition_counts = np.zeros((2, 2, 2, 2), dtype=float)

    # initial state, final state, state
    dwell_times = np.zeros((2, 2, 2), dtype=float)

    # get the path samples and track some stats
    T = 1.0
    npaths = 1000
    for i in range(npaths):
        t_counts = np.zeros((2, 2), dtype=float)
        d_times = np.zeros(2, dtype=float)
        initial_state = random_category(initial_distn)
        t = 0.0
        state = initial_state
        while t < T and state in (0, 1):
            rate = -Q[state, state]
            dwell = random.expovariate(rate)
            t_next = t + dwell
            if t_next > T:
                d_times[state] += T - t
            else:
                s_next = random_category(P[state])
                if s_next in (0, 1):
                    d_times[state] += dwell
                    t_counts[state, s_next] += 1
                    state = s_next
            t = t_next
        final_state = state
        if final_state not in (0, 1):
            continue
        transition_counts[initial_state, final_state] += t_counts
        dwell_times[initial_state, final_state] += d_times


    # report stuff
    print >> out, 'rate a:', a
    print >> out, 'rate b:', b
    print >> out, 'rate c:', c
    print >> out
    print >> out, 'U:'
    print >> out, U
    print >> out
    print >> out, 'w:',
    print >> out, w
    print >> out
    print >> out, 'V:'
    print >> out, V
    print >> out
    print >> out, 'U V (should be the 3x3 identity):'
    print >> out, np.dot(U, V)
    print >> out
    print >> out, 'Q (the rate matrix):'
    print >> out, Q
    print >> out
    print >> out, 'U diag(w) V (should be the same as Q):'
    print >> out, np.dot(U * w, V)
    print >> out
    print >> out
    print >> out, '=== path sample statistics ==='
    print >> out
    for a in (0, 1):
        for b in (0, 1):
            print >> out, 'path', a, '--->', b
            print >> out, 'transition counts:'
            print >> out, transition_counts[a, b]
            print >> out, 'dwell times:'
            print >> out, dwell_times[a, b]
            print >> out

    # show the result
    return out.getvalue()

