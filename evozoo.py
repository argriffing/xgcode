"""
The intention is to collect a bunch of parameterized rate matrices.

A focus is on site-independent mutation and site-dependent selection,
using a specific formula to obtain the mutation-selection balance
rate matrix given the mutation rate matrix and the
per-state selection parameters.
"""

import numpy as np

import mrate


# The df member variable is the number of degrees of freedom.
# The X parameter for get_rate_matrix is an ndarray which
# should have a number of selection 'energy' values
# that is equal to the degrees of freedom df.
# The idea is that because energies can take any float value,
# it is straightforward to use a numerical search to minimize over
# such an unconstrained parameter space.
# For general selection the number of such unconstrained energy parameters
# is one fewer than the number of states of the process,
# but for some special selection constraints
# this number of degrees of freedom is reduced.

# For most of these processes,
# the mutation process will have each rate the same and equal to 1.

def energies_to_distn(v):
    p = np.exp(-v)
    return p / np.sum(p)

def get_uniform_distn(nstates):
    return np.ones(nstates, dtype=float) / nstates

class Process:
    def _check_params(self, X):
        if len(X) != self.get_df():
            raise ValueError('mismatch of number of degrees of freedom')

# In the following sequence of processes,
# the _x_y_z describe
# the number of states per site,
# the number of sites,
# and the number of selection parameters.

class Hypercube_2_3_0(Process):
    """
    This is a cube.
    """
    def get_df(self):
        return 0
    def get_nstates(self):
        return 8
    def get_distn(self, X):
        self._check_params(X)
        return get_uniform_distn(self.get_nstates())
    def get_rate_matrix(self, X):
        self._check_params(X)
        return np.array([
            [-3, 1, 1, 0, 1, 0, 0, 0],
            [1, -3, 0, 1, 0, 1, 0, 0],
            [1, 0, -3, 1, 0, 0, 1, 0],
            [0, 1, 1, -3, 0, 0, 0, 1],
            [1, 0, 0, 0, -3, 1, 1, 0],
            [0, 1, 0, 0, 1, -3, 0, 1],
            [0, 0, 1, 0, 1, 0, -3, 1],
            [0, 0, 0, 1, 0, 1, 1, -3]], dtype=float)

class Coil_2_3_0(Process):
    """
    This is a maximal induced cycle in a cube.
    In graph theory jargon it is also known as coil-in-the-box.
    """
    def get_df(self):
        return 0
    def get_nstates(self):
        return 6
    def get_distn(self, X):
        self._check_params(X)
        return get_uniform_distn(self.get_nstates())
    def get_rate_matrix(self, X):
        self._check_params(X)
        return np.array([
            [-2, 1, 0, 0, 0, 1],
            [1, -2, 1, 0, 0, 0],
            [0, 1, -2, 1, 0, 0],
            [0, 0, 1, -2, 1, 0],
            [0, 0, 0, 1, -2, 1],
            [1, 0, 0, 0, 1, -2]], dtype=float)

class Snake_2_3_0(Process):
    """
    This is a maximal induced path in a cube.
    In graph theory jargon it is also known as snake-in-the-box.
    """
    def get_df(self):
        return 0
    def get_nstates(self):
        return 5
    def get_distn(self, X):
        self._check_params(X)
        return get_uniform_distn(self.get_nstates())
    def get_rate_matrix(self, X):
        self._check_params(X)
        return np.array([
            [-1, 1, 0, 0, 0],
            [1, -2, 1, 0, 0],
            [0, 1, -2, 1, 0],
            [0, 0, 1, -2, 1],
            [0, 0, 0, 1, -1]], dtype=float)


# The following sequence of processes each have a single selection parameter.
# This parameter is unconstrained and controls the log ratio of the
# probabiliies of alternating states, for a specific meaning of alternating.
# For the cycle and path, alternating states mean what you would expect.
# For the hypercube, alternating states mean states of even vs. odd parity
# of the binary strings corresponding to the states,
# where adjacent states have hamming distance 1 between their
# corresponding binary strings.
# Note that all of these meanings of alternating give induced
# graphs that become nearly disconnected when the probability ratio is
# far from 1, in the sense that each path between high-probability
# states must pass through a low-probability state.


class AlternatingHypercube_2_3_1(Process):
    def get_df(self):
        return 1
    def get_nstates(self):
        return 8
    def get_distn(self, X):
        self._check_params(X)
        return energies_to_distn(self._get_energies(X))
    def _get_energies(self, X):
        self._check_params(X)
        g, = X
        return np.array([g, 0, 0, g, 0, g, g, 0])
    def get_rate_matrix(self, X):
        self._check_params(X)
        M = np.array([
            [-3, 1, 1, 0, 1, 0, 0, 0],
            [1, -3, 0, 1, 0, 1, 0, 0],
            [1, 0, -3, 1, 0, 0, 1, 0],
            [0, 1, 1, -3, 0, 0, 0, 1],
            [1, 0, 0, 0, -3, 1, 1, 0],
            [0, 1, 0, 0, 1, -3, 0, 1],
            [0, 0, 1, 0, 1, 0, -3, 1],
            [0, 0, 0, 1, 0, 1, 1, -3]], dtype=float)
        u = np.zeros(self.get_nstates())
        v = self._get_energies(X)
        return mrate.to_gtr_hb_known_energies(M, u, v)

class AlternatingCoil_2_3_1(Process):
    def get_df(self):
        return 1
    def get_nstates(self):
        return 6
    def get_distn(self, X):
        self._check_params(X)
        return energies_to_distn(self._get_energies(X))
    def _get_energies(self, X):
        self._check_params(X)
        g, = X
        return np.array([g, 0, g, 0, g, 0])
    def get_rate_matrix(self, X):
        self._check_params(X)
        M =  np.array([
            [-2, 1, 0, 0, 0, 1],
            [1, -2, 1, 0, 0, 0],
            [0, 1, -2, 1, 0, 0],
            [0, 0, 1, -2, 1, 0],
            [0, 0, 0, 1, -2, 1],
            [1, 0, 0, 0, 1, -2]], dtype=float)
        u = np.zeros(self.get_nstates())
        v = self._get_energies(X)
        return mrate.to_gtr_hb_known_energies(M, u, v)

class AlternatingSnake_2_3_1(Process):
    def get_df(self):
        return 1
    def get_nstates(self):
        return 5
    def get_distn(self, X):
        self._check_params(X)
        return energies_to_distn(self._get_energies(X))
    def _get_energies(self, X):
        self._check_params(X)
        g, = X
        return np.array([g, 0, g, 0, g])
    def get_rate_matrix(self, X):
        self._check_params(X)
        M = np.array([
            [-1, 1, 0, 0, 0],
            [1, -2, 1, 0, 0],
            [0, 1, -2, 1, 0],
            [0, 0, 1, -2, 1],
            [0, 0, 0, 1, -1]], dtype=float)
        u = np.zeros(self.get_nstates())
        v = self._get_energies(X)
        return mrate.to_gtr_hb_known_energies(M, u, v)

