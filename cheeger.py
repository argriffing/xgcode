"""
Do things related to isoperimetric constants of Markov processes.

The Markov processes are assumed to be
continuous-time, irreducible, reversible, finite-state.
"""

from StringIO import StringIO
import itertools
from itertools import product

import numpy as np

import iterutils


def _get_cheeger_constant(R, v):
    """
    This is also known as the second isoperimetric constant.
    @param R: a reversible rate matrix
    @param v: stationary distribution
    @return: the second isoperimetric constant
    """
    n = len(v)
    I2 = None
    for A_tuple in iterutils.powerset(range(n)):
        # define the vertex set and its complement
        A = set(A_tuple)
        B = set(range(n)) - A
        A_measure = sum(v[i] for i in A)
        B_measure = sum(v[i] for i in B)
        if A_measure and B_measure:
            boundary_measure = sum(v[i]*R[i,j] for i, j in product(A, B))
            A_connectivity = boundary_measure / A_measure
            B_connectivity = boundary_measure / B_measure
            connectivity = max(A_connectivity, B_connectivity)
            if I2 is None or connectivity < I2:
                I2 = connectivity
    return I2

def get_cheeger_bounds(R, v):
    """
    @param R: a reversible rate matrix
    @param v: stationary distribution
    @return: low, cheeger, high
    """
    cheeger = _get_cheeger_constant(R, v)
    low = 0.5 * (cheeger**2) / -min(np.diag(R))
    high = 2 * cheeger
    return low, cheeger, high

