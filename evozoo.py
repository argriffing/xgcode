"""
The intention is to collect a bunch of parameterized rate matrices.

A focus is on site-independent mutation and site-dependent selection,
using a specific formula to obtain the mutation-selection balance
rate matrix given the mutation rate matrix and the
per-state selection parameters.
"""

import numpy as np
import gmpy

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

def coil_in_the_box(d):
    """
    OEIS A000937
    @param d: hypercube dimension
    @return: edges in the longest inducible cycle
    """
    d_to_edge_count = {
            2 : 4,
            3 : 6,
            4 : 8,
            5 : 14,
            6 : 26,
            7 : 48,
            }
    dstr = '(d = %d)' % d
    if d < 2:
        raise ValueError('coil-in-the-box does not exist ' + dstr)
    if d not in d_to_edge_count:
        raise ValueError('coil-in-the-box has not been solved ' + dstr)
    return d_to_edge_count[d]

def snake_in_the_box(d):
    """
    OEIS A099155
    @param d: hypercube dimension
    @return: edges in the longest inducible path
    """
    d_to_edge_count = {
            1 : 1,
            2 : 2,
            3 : 4,
            4 : 7,
            5 : 13,
            6 : 26,
            7 : 50,
            }
    dstr = '(d = %d)' % d
    if d < 2:
        raise ValueError('snake-in-the-box does not exist ' + dstr)
    if d not in d_to_edge_count:
        raise ValueError('snake-in-the-box has not been solved ' + dstr)
    return d_to_edge_count[d]

def energies_to_distn(v):
    p = np.exp(-v)
    return p / np.sum(p)

def get_uniform_distn(nstates):
    return np.ones(nstates, dtype=float) / nstates

class Process:
    def _check_params(self, X):
        if len(X) != self.get_df():
            raise ValueError('mismatch of number of degrees of freedom')

# The following cycle and path have no selection parameters.
# The N refers to the number of states.
# The 0 refers to the degrees of freedom of selection.

class Cycle_N_0(Process):
    def __init__(self, nstates):
        self.nstates = nstates
    def get_df(self):
        return 0
    def get_nstates(self):
        return self.nstates
    def get_distn(self, X=None):
        if X is not None:
            self._check_params(X)
        return get_uniform_distn(self.get_nstates())
    def get_rate_matrix(self, X=None):
        if X is not None:
            self._check_params(X)
        n = self.get_nstates()
        Q = np.zeros((n, n))
        for i in range(n):
            Q[i,i] = -2
        for i in range(n-1):
            Q[i,i+1] = 1
            Q[i+1,i] = 1
        Q[0,n-1] = 1
        Q[n-1,0] = 1
        return Q

class Path_N_0(Process):
    def __init__(self, nstates):
        self.nstates = nstates
    def get_df(self):
        return 0
    def get_nstates(self):
        return self.nstates
    def get_distn(self, X=None):
        if X is not None:
            self._check_params(X)
        return get_uniform_distn(self.get_nstates())
    def get_rate_matrix(self, X=None):
        if X is not None:
            self._check_params(X)
        n = self.get_nstates()
        Q = np.zeros((n, n))
        for i in range(n):
            Q[i,i] = -2
        for i in range(n-1):
            Q[i,i+1] = 1
            Q[i+1,i] = 1
        Q[0,0] = -1
        Q[n-1,n-1] = -1
        return Q

# This is a hypercube graph with 2^d states.
# It has no selection parameters.
# Because it is a hypercube, the number of states per
# dimension is hardcoded to two.

class Hypercube_d_0(Process):
    def __init__(self, d):
        self.d = d
    def get_df(self):
        return 0
    def get_nstates(self):
        return 2**self.d
    def get_distn(self, X=None):
        if X is not None:
            self._check_params(X)
        return get_uniform_distn(self.get_nstates())
    def get_rate_matrix(self, X=None):
        if X is not None:
            self._check_params(X)
        d = self.d
        n = self.get_nstates()
        Q = np.zeros((n, n))
        for a in range(n):
            for b in range(n):
                dist = gmpy.hamdist(a, b)
                if dist == 0:
                    Q[a, b] = -d
                elif dist == 1:
                    Q[a, b] = 1
        return Q


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


# The following are single selection parameter processes.

class _Alternating(Process):
    def get_df(self):
        return 1
    def get_distn(self, X):
        return energies_to_distn(self.get_energies(X))

class AlternatingHypercube_d_1(_Alternating):
    def __init__(self, d):
        self.d = d
    def get_nstates(self):
        return 2**self.d
    def get_energies(self, X):
        self._check_params(X)
        g, = X
        popcounts = [gmpy.popcount(i) for i in range(self.get_nstates())]
        return np.array([g if p%2 else 0 for p in popcounts])
    def get_rate_matrix(self, X):
        self._check_params(X)
        M = Hypercube_d_0(self.d).get_rate_matrix()
        u = np.zeros(self.get_nstates())
        v = self.get_energies(X)
        return mrate.to_gtr_hb_known_energies(M, u, v)

class AlternatingCycle_N_1(_Alternating):
    def __init__(self, nstates):
        self.nstates = nstates
    def get_nstates(self):
        return self.nstates
    def get_energies(self, X):
        self._check_params(X)
        g, = X
        return np.array([g if i%2 else 0 for i in range(self.get_nstates())])
    def get_rate_matrix(self, X):
        self._check_params(X)
        n = self.get_nstates()
        M = Cycle_N_0(n).get_rate_matrix()
        u = np.zeros(n)
        v = self.get_energies(X)
        return mrate.to_gtr_hb_known_energies(M, u, v)

class AlternatingPath_N_1(_Alternating):
    def __init__(self, nstates):
        self.nstates = nstates
    def get_nstates(self):
        return self.nstates
    def get_energies(self, X):
        self._check_params(X)
        g, = X
        return np.array([g if i%2 else 0 for i in range(self.get_nstates())])
    def get_rate_matrix(self, X):
        self._check_params(X)
        n = self.get_nstates()
        M = Path_N_0(n).get_rate_matrix()
        u = np.zeros(n)
        v = self.get_energies(X)
        return mrate.to_gtr_hb_known_energies(M, u, v)


# The following are fully parameterized selection processes.

class _General(Process):
    def get_df(self):
        return self.get_nstates() - 1
    def get_distn(self, X):
        return energies_to_distn(self.get_energies(X))
    def get_energies(self, X):
        self._check_params(X)
        return np.hstack(([0], X))

class GeneralHypercube_d_full(_General):
    def __init__(self, d):
        self.d = d
    def get_nstates(self):
        return 2**self.d
    def get_rate_matrix(self, X):
        self._check_params(X)
        M = Hypercube_d_0(self.d).get_rate_matrix()
        u = np.zeros(self.get_nstates())
        v = self.get_energies(X)
        return mrate.to_gtr_hb_known_energies(M, u, v)

class GeneralCycle_N_full(_General):
    def __init__(self, nstates):
        self.nstates = nstates
    def get_nstates(self):
        return self.nstates
    def get_rate_matrix(self, X):
        self._check_params(X)
        n = self.get_nstates()
        M = Cycle_N_0(n).get_rate_matrix()
        u = np.zeros(n)
        v = self.get_energies(X)
        return mrate.to_gtr_hb_known_energies(M, u, v)

class GeneralPath_N_full(_General):
    def __init__(self, nstates):
        self.nstates = nstates
    def get_nstates(self):
        return self.nstates
    def get_rate_matrix(self, X):
        self._check_params(X)
        n = self.get_nstates()
        M = Path_N_0(n).get_rate_matrix()
        u = np.zeros(n)
        v = self.get_energies(X)
        return mrate.to_gtr_hb_known_energies(M, u, v)


# These are coil-in-the-box and snake-in-the-box.
# http://en.wikipedia.org/wiki/Snake-in-the-box
# No selection parameters.
# The d is the number of dimensions of the embedding hypercube.

class Coil_d_0(Cycle_N_0):
    def __init__(self, d):
        nstates = coil_in_the_box(d)
        Cycle_N_0.__init__(self, nstates)

class Snake_d_0(Path_N_0):
    def __init__(self, d):
        nstates = snake_in_the_box(d) + 1
        Path_N_0.__init__(self, nstates)


# These are coil-in-the-box and snake-in-the-box with a selection parameter.
# The d is the number of dimensions of the embedding hypercube.

class AlternatingCoil_d_1(AlternatingCycle_N_1):
    def __init__(self, d):
        nstates = coil_in_the_box(d)
        AlternatingCycle_N_1.__init__(self, nstates)

class AlternatingSnake_d_1(AlternatingPath_N_1):
    def __init__(self, d):
        nstates = snake_in_the_box(d) + 1
        AlternatingPath_N_1.__init__(self, nstates)


# Fully parameterized.

class GeneralCoil_d_full(GeneralCycle_N_full):
    def __init__(self, d):
        nstates = coil_in_the_box(d)
        GeneralCycle_N_full.__init__(self, nstates)

class GeneralSnake_d_full(GeneralPath_N_full):
    def __init__(self, d):
        nstates = snake_in_the_box(d) + 1
        GeneralPath_N_full.__init__(self, nstates)


# A one-off.
# A single pair of corners of the hypercube is distinguished.
# When d=3 this interpolates between the 3-cube and its largest induced cycle.

class DistinguishedCornerPairHypercube_d_1(Process):
    def __init__(self, d):
        self.d = d
    def get_df(self):
        return 1
    def get_distn(self, X):
        return energies_to_distn(self.get_energies(X))
    def get_nstates(self):
        return 2 ** self.d
    def get_energies(self, X):
        self._check_params(X)
        g, = X
        popcounts = [gmpy.popcount(i) for i in range(self.get_nstates())]
        return np.array([g if p in (0, self.d) else 0 for p in popcounts])
    def get_rate_matrix(self, X):
        self._check_params(X)
        M = Hypercube_d_0(self.d).get_rate_matrix()
        u = np.zeros(self.get_nstates())
        v = self.get_energies(X)
        return mrate.to_gtr_hb_known_energies(M, u, v)

