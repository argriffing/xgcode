"""
Utility functions for finite Markov chains.
"""

import unittest

import numpy as np
from scipy import linalg

import MatrixUtil


def inverse_permutation(p):
    p_inv = [0] * len(p)
    for i, j in enumerate(p):
        p_inv[j] = i
    return np.array(p_inv)

def get_conditional_transition_matrix(P, plain, forbid, target):
    """
    Get a conditional transition matrix.
    The new transition matrix is conditioned on reaching a state
    in the sequence of target states before reaching a state
    in the sequence of forbidden states.
    @param P: transition matrix
    @param plain: sequence of plain state indices
    @param forbid: sequence of forbidden state indices
    @param target: sequence of target state indices
    @return: a conditional transition matrix conformant with P
    """
    # check that P is really a transition matrix
    MatrixUtil.assert_transition_matrix(P)
    # define some state lists
    special = np.hstack((forbid, target))
    states = np.hstack((plain, special))
    # check that the index sequences match the size of P
    if sorted(states) != range(len(P)):
        raise ValueError('P is not conformant with the index sequences')
    # Q is the part of the transition matrix that gives
    # transition probabilities from transient to transient states.
    # In general it is not row stochastic.
    Q = P[plain, :][:, plain]
    # R is the part of the transition matrix that gives
    # transition probabilities from transient to absorbing states.
    # In general it is not row stochastic.
    R = P[plain, :][:, special]
    # Note that Q and R completely define an absorbing process,
    # because there are no transitions away from absorbing states.
    # Each row of B is a probability distribution.
    B = linalg.solve(np.eye(len(plain)) - Q, R)
    # Use the distributions in B to compute weights
    # to apply to the transitions.
    c = np.hstack((np.zeros_like(forbid), np.ones_like(target)))
    w = np.hstack((np.dot(B, c), c))
    # Make a new transition matrix that is conditioned
    # on hitting a target state before hitting a forbidden state.
    # Permute the weights to conform to the original matrix indices.
    # Rescale rows to restore row-stochasticity.
    H = P * w[inverse_permutation(states)]
    v = np.sum(H, axis=1)
    H /= v[:, np.newaxis]
    return H

def get_absorption_time(P, plain, absorbing):
    """
    Get expected times to absorption.
    Note that if an index is indicated as absorbing by its presence
    in the sequence of absorbing state indices,
    then it will be treated as absorbing
    even if the transition matrix P indicates otherwise.
    @param P: transition matrix
    @param plain: sequence of plain state indices
    @param absorbing: sequence of absorbing state indices
    @return: expected times to absorption or 0 from absorbing states
    """
    # check that P is really a transition matrix
    MatrixUtil.assert_transition_matrix(P)
    # define some state lists
    states = np.hstack((plain, absorbing))
    # check that the index sequences match the size of P
    if sorted(states) != range(len(P)):
        raise ValueError('P is not conformant with the index sequences')
    # compute the time to absorption
    Q = P[plain, :][:, plain]
    c = np.ones_like(plain)
    tplain = linalg.solve(np.eye(len(plain)) - Q, c)
    t = np.hstack((tplain, np.zeros_like(absorbing)))
    return t[inverse_permutation(states)]

def get_absorption_variance(P, plain, absorbing):
    """
    Get expected times to absorption.
    Note that if an index is indicated as absorbing by its presence
    in the sequence of absorbing state indices,
    then it will be treated as absorbing
    even if the transition matrix P indicates otherwise.
    @param P: transition matrix
    @param plain: sequence of plain state indices
    @param absorbing: sequence of absorbing state indices
    @return: variance of times to absorption or 0 from absorbing states
    """
    # check that P is really a transition matrix
    MatrixUtil.assert_transition_matrix(P)
    # define some state lists
    states = np.hstack((plain, absorbing))
    # check that the index sequences match the size of P
    if sorted(states) != range(len(P)):
        raise ValueError('P is not conformant with the index sequences')
    # compute the time to absorption
    Q = P[plain, :][:, plain]
    c = np.ones_like(plain)
    I = np.eye(len(plain))
    t = linalg.solve(I - Q, c)
    # compute the variance
    vplain = 2*linalg.solve(I - Q, t) - t*(t+1)
    v = np.hstack((vplain, np.zeros_like(absorbing)))
    return v[inverse_permutation(states)]

class TestHittingTime(unittest.TestCase):
    def test_hitting_time(self):
        pass
    def test_absorption_time(self):
        pass

if __name__ == '__main__':
    unittest.main()

