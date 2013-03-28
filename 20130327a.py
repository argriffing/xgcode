"""
Check properties of an ising-like chain model x0--x1--x2.
"""

from StringIO import StringIO
import itertools
from collections import defaultdict

import numpy as np

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

def get_energy(state):
    """
    This computes a hamiltonian function.
    It is not a problem if this function is slow.
    @param state: a configuration
    @return: the energy of the state
    """
    v = np.array(state, dtype=float) - 0.5
    zeroth_order = np.sum(v)
    first_order = np.sum(v[1:] * v[:-1])
    return zeroth_order + first_order

def is_wildmatch(key, state):
    for k, s in zip(key, state):
        if k != '?' and k != s:
            return False
    return True

def keystr(key):
    return '(%s, %s, %s)' % key

def get_response_content(fs):

    # set up print options
    np.set_printoptions(
            linewidth=1000000,
            threshold=1000000,
            )
    out = StringIO()

    # Define the states.
    states = list(itertools.product((0, 1), repeat=3))

    # Define the joint distribution.
    energies = np.array([get_energy(state) for state in states])
    unnormalized_distn = np.exp(-energies)
    distn = unnormalized_distn / np.sum(unnormalized_distn)
    print >> out, 'energies of (x0, x1, x2):'
    for state, energy in zip(states, energies):
        print >> out, state, energy
    print >> out
    print >> out, 'distribution of (x0, x1, x2):'
    for state, p in zip(states, distn):
        print >> out, state, p
    print >> out

    # Define some marginal probabilities.
    marginals = defaultdict(float)
    for key in itertools.product((0, 1, '?'), repeat=3):
        for state, p in zip(states, distn):
            if is_wildmatch(key, state):
                marginals[key] += p

    # Print the marginals.
    print >> out, 'marginal probabilities:'
    for key, p in sorted(marginals.items()):
        print >> out, keystr(key), p
    print >> out

    # Print a distribution that assumes a Markov structure.
    print >> out, 'P(x0) * P(x1 | x0) * P(x2 | x1) :'
    for state in states:
        x0, x1, x2 = state
        p = 1.0

        # P(x0)
        p *= marginals[(x0, '?', '?')]

        # P(x1 | x0)
        p *= marginals[(x0, x1, '?')]
        p /= marginals[(x0, '?', '?')]

        # P(x2 | x1)
        p *= marginals[('?', x1, x2)]
        p /= marginals[('?', x1, '?')]

        print >> out, state, p
    print >> out

    # Print a distribution that uses a truncated factorization.
    print >> out, 'truncated probability decomposition:'
    for state in states:
        x0, x1, x2 = state
        xx = '?'
        p = 1.0

        #  P(x0) * P(x1) * P(x2)
        p *= marginals[(x0, xx, xx)]
        p *= marginals[(xx, x1, xx)]
        p *= marginals[(xx, xx, x2)]

        # P(x0, x1) / ( P(x0) * P(x1) )
        p *= marginals[(x0, x1, xx)]
        p /= marginals[(x0, xx, xx)]
        p /= marginals[(xx, x1, xx)]

        # P(x1, x2) / ( P(x1) * P(x2) )
        p *= marginals[(xx, x1, x2)]
        p /= marginals[(xx, x1, xx)]
        p /= marginals[(xx, xx, x2)]

        print >> out, state, p
    print >> out
    
    """
    #
    # compute some conditional distributions
    #
    x0_distn_given_x1_n1 = np.zeros(2)
    for state, p in zip(states, distn):
        if state[1] == -1:
            if state[0] == -1:
                x0_distn_given_x1_n1[0] += p
            elif state[0] == 1:
                x0_distn_given_x1_n1[1] += p
    x0_distn_given_x1_n1 /= np.sum(x0_distn_given_x1_n1)
    #
    x2_distn_given_x1_n1 = np.zeros(2)
    for state, p in zip(states, distn):
        if state[1] == -1:
            if state[2] == -1:
                x2_distn_given_x1_n1[0] += p
            elif state[2] == 1:
                x2_distn_given_x1_n1[1] += p
    x2_distn_given_x1_n1 /= np.sum(x2_distn_given_x1_n1)
    #
    x0_distn_given_x1_p1 = np.zeros(2)
    for state, p in zip(states, distn):
        if state[1] == 1:
            if state[0] == -1:
                x0_distn_given_x1_p1[0] += p
            elif state[0] == 1:
                x0_distn_given_x1_p1[1] += p
    x0_distn_given_x1_p1 /= np.sum(x0_distn_given_x1_p1)
    #
    x2_distn_given_x1_p1 = np.zeros(2)
    for state, p in zip(states, distn):
        if state[1] == 1:
            if state[2] == -1:
                x2_distn_given_x1_p1[0] += p
            elif state[2] == 1:
                x2_distn_given_x1_p1[1] += p
    x2_distn_given_x1_p1 /= np.sum(x2_distn_given_x1_p1)
    #
    # compute a couple more conditional distributions
    #
    """

    # show the result
    return out.getvalue()

