"""
Look at the algebraic connectivity of random subgraphs of hypercubes.

The graphs are unweighted.
"""

from StringIO import StringIO
import math
import itertools
import random

import numpy as np
import scipy
from scipy import linalg

import Form
import FormOut
from MatrixUtil import ndot
import mrate
import combobreaker
import iterutils

def get_form():
    """
    @return: the body of a form
    """
    # define the form objects
    form_objects = [
            Form.Integer('ndim', 'hypercube dimension',
                5, low=2, high=8),
            Form.Float('pkeep', 'probability to keep each vertex',
                0.5, low_exclusive=0, high_inclusive=1)]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_local_cheeger_ratio_bound(Q):
    """
    Return an upper bound on the local Cheeger ratio.
    The local Cheeger ratio is an upper bound on the Cheeger ratio.
    The Cheeger ratio is half of the upper bound on the algebraic connectivity.
    This function is randomly jittered to avoid symmetry.
    @param Q: a rate matrix
    @return: an upper bound on the cheeger ratio
    """
    jitter = 0.001
    jitter_low = jitter - jitter / 2
    jitter_high = jitter + jitter / 2
    n = len(Q)
    A_jittered = np.copy(Q) - np.diag(np.diag(Q))
    for i in range(n):
        for j in range(n):
            if i < j:
                r = math.exp(random.uniform(jitter_low, jitter_high))
                A_jittered[i, j] *= r
                A_jittered[j, i] *= r
    Q_jittered = A_jittered - np.diag(np.sum(A_jittered, axis=1))
    W, V = scipy.linalg.eigh(
            Q_jittered,
            eigvals=(n-2, n-1))
    fiedler = V.T[0]
    neg_values = [x for x in fiedler if x < 0]
    pos_values = [x for x in fiedler if x >= 0]
    vertex_size = min(len(neg_values), len(pos_values))
    boundary_size = 0
    for i in range(n):
        for j in range(n):
            if i < j:
                if Q[i,j] and (fiedler[i] * fiedler[j] < 0):
                    boundary_size += 1
    ratio_bound = float(boundary_size) / float(vertex_size)
    return ratio_bound

def get_real_cheeger(Q):
    """
    Get the real cheeger constant.
    This is not a bound or an approximation.
    And it is NP hard to compute.
    """
    n = len(Q)
    all_vertices = set(range(n))
    min_ratio = None
    for vertices in iterutils.powerset(range(n)):
        if not (0 < len(vertices) <= n/2):
            continue
        vset = set(vertices)
        complement = all_vertices - vset
        volume = len(vset)
        boundary_size = 0
        for v in vset:
            for x in complement:
                if Q[v, x]:
                    boundary_size += 1
        ratio = boundary_size / float(volume)
        if (min_ratio is None) or (ratio < min_ratio):
            min_ratio = ratio
    return min_ratio

def get_detailed_balance_error(Q):
    """
    @param Q: a rate matrix
    @return: a number that should be near zero if detailed balance is satisfied
    """
    p = mrate.R_to_distn(Q)
    errors = []
    nstates = len(Q)
    for i in range(nstates):
        for j in range(nstates):
            error = p[i] * Q[i, j] - p[j] * Q[j, i]
            errors.append(error)
    return max(abs(x) for x in errors)

def get_rate_matrix_summary(Q):
    out = StringIO()
    Q_t = mrate.R_to_relaxation_time(Q)
    Q_cheeger_bound = get_local_cheeger_ratio_bound(Q)
    if len(Q) < 16:
        real_cheeger = get_real_cheeger(Q)
        cheeger_string = str(real_cheeger)
    else:
        cheeger_string = 'takes too long to compute'
    print >> out, 'rate matrix:'
    print >> out, Q
    print >> out
    print >> out, 'algebraic connectivity (all hypothesized to be <= 2)'
    print >> out, 1 / Q_t
    print >> out
    print >> out, 'local cheeger ratio bound (all hypothesized to be <= 1):'
    print >> out, Q_cheeger_bound
    print >> out
    print >> out, 'actual Cheeger constant (all hypothesized to be <= 1):'
    print >> out, cheeger_string
    print >> out
    return out.getvalue().rstrip()

def get_unweighted_path(nstates):
    """
    @return: a rate matrix
    """
    A = np.zeros((nstates, nstates))
    for i in range(nstates-1):
        a, b = i, (i+1) % nstates
        A[a,b] = 1
        A[b,a] = 1
    Q = A - np.diag(np.sum(A, axis=1))
    return Q

def get_hamming_graph_adjacency(nresidues, nsites):
    """
    The number of states is (nresidues ** nsites).
    @return: an adjacency matrix
    """
    nstates = nresidues**nsites
    A = np.zeros((nstates, nstates))
    for alpha in itertools.product(range(nresidues), repeat=nsites):
        for beta in itertools.product(range(nresidues), repeat=nsites):
            alpha_index = sum(alpha[i]*(nresidues ** i) for i in range(nsites))
            beta_index = sum(beta[i]*(nresidues ** i) for i in range(nsites))
            hamming_dist = sum(1 for a, b in zip(alpha, beta) if a != b)
            if hamming_dist == 1:
                A[alpha_index, beta_index] = 1
    return A

def get_hypercube_adjacency(d):
    """
    The number of states is (2 ** d).
    @param d: the number of dimensions
    """
    return get_hamming_graph_adjacency(2, d)

def get_largest_component(A):
    """
    Return the adjacency matrix of the largest connected component.
    @param A: a connected unweighted adjacency matrix
    @return: None a connected unweighted adjacency matrix 
    """
    n = len(A)
    A_sparse = [[x for x in range(n) if A[v,x]] for v in range(n)]
    visited = set()
    component_vertex_sets = []
    for seed in range(n):
        if seed in visited:
            continue
        component_vertex_set = set()
        shell = {seed}
        visited.add(seed)
        while shell:
            next_shell = set()
            for v in shell:
                component_vertex_set.add(v)
                for x in A_sparse[v]:
                    if x not in visited:
                        next_shell.add(x)
                        visited.add(x)
            shell = next_shell
        component_vertex_sets.append(component_vertex_set)
    pairs = [(len(x), x) for x in component_vertex_sets]
    best_size, best_set = max(pairs)
    if not best_size:
        return None
    big_to_small = dict((v, i) for i, v in enumerate(sorted(best_set)))
    A_component = np.zeros((best_size, best_size))
    for v_big in best_set:
        v_small = big_to_small[v_big]
        for x_big in A_sparse[v_big]:
            x_small = big_to_small[x_big]
            A_component[v_small, x_small] = 1
    return A_component

class MyCallback:
    def __init__(self, ndim, pkeep):
        self.ndim = ndim
        self.pkeep = pkeep
        self.connectivities = []
        self.cheegers = []
        self.max_connectivity = None
        self.max_connectivity_graph = None
        self.min_connectivity = None
        self.min_connectivity_graph = None
        self.max_cheeger = None
        self.max_cheeger_graph = None
        self.min_cheeger = None
        self.min_cheeger_graph = None
    def __call__(self):
        nstates = 2 ** self.ndim
        # get a hypercube with the given number of dimensions
        A_full = get_hypercube_adjacency(self.ndim)
        # sample a random hypercube by removing a subset of vertices
        keep_list = [v for v in range(nstates) if random.random() < self.pkeep]
        nkeep = len(keep_list)
        if nkeep < 2:
            return
        A_reduced = np.zeros((nkeep, nkeep))
        for i, v in enumerate(keep_list):
            for j, x in enumerate(keep_list):
                A_reduced[i, j] = A_full[v, x]
        # get the largest connected component of the random induced subgraph
        A_component = get_largest_component(A_reduced)
        if A_component is None:
            return
        ncomponent = len(A_component)
        if ncomponent < 2:
            return
        Q_component = A_component - np.diag(np.sum(A_component, axis=1))
        w = scipy.linalg.eigvalsh(
                Q_component,
                eigvals=(ncomponent-2, ncomponent-1))
        neg_rate = w[0]
        # compute the connectivity and the cheeger bound
        connectivity = -neg_rate
        cheeger = get_local_cheeger_ratio_bound(Q_component)
        # update the induced subgraphs with the extreme statistics
        self._update_min_connectivity(Q_component, connectivity)
        self._update_max_connectivity(Q_component, connectivity)
        self._update_min_cheeger(Q_component, cheeger)
        self._update_max_cheeger(Q_component, cheeger)
        # add the connectivity and the cheeger bounds
        self.connectivities.append(connectivity)
        self.cheegers.append(cheeger)
    def _update_min_connectivity(self, Q_component, mu):
        if (not self.connectivities) or (self.min_connectivity > mu):
            self.min_connectivity = mu
            self.min_connectivity_graph = Q_component
    def _update_max_connectivity(self, Q_component, mu):
        if (not self.connectivities) or (self.max_connectivity < mu):
            self.max_connectivity = mu
            self.max_connectivity_graph = Q_component
    def _update_min_cheeger(self, Q_component, cheeger):
        if (not self.cheegers) or (self.min_cheeger > cheeger):
            self.min_cheeger = cheeger
            self.min_cheeger_graph = Q_component
    def _update_max_cheeger(self, Q_component, cheeger):
        if (not self.cheegers) or (self.max_cheeger < cheeger):
            self.max_cheeger = cheeger
            self.max_cheeger_graph = Q_component
    def __str__(self):
        snake_lengths = [0, 1, 2, 4, 7, 13, 26, 50, 96, 188, 348, 640, 1238]
        nsnake = snake_lengths[self.ndim] + 1
        Q_snake = get_unweighted_path(nsnake)
        snake_connect = -scipy.linalg.eigvalsh(
                Q_snake,
                eigvals=(nsnake-2, nsnake-1))[0]
        A_full = get_hypercube_adjacency(self.ndim)
        Q_full = A_full - np.diag(np.sum(A_full, axis=1))
        out = StringIO()
        print >> out
        print >> out, 'full %d-cube:' % self.ndim
        print >> out
        print >> out, get_rate_matrix_summary(Q_full)
        print >> out
        print >> out
        print >> out, 'min algebraic connectivity induced subgraph:'
        print >> out
        print >> out, get_rate_matrix_summary(self.min_connectivity_graph)
        print >> out
        print >> out
        print >> out, 'max algebraic connectivity induced subgraph:'
        print >> out
        print >> out, get_rate_matrix_summary(self.max_connectivity_graph)
        print >> out
        print >> out
        print >> out, 'min local cheeger ratio bound induced subgraph:'
        print >> out
        print >> out, get_rate_matrix_summary(self.min_cheeger_graph)
        print >> out
        print >> out
        print >> out, 'max local cheeger ratio bound induced subgraph:'
        print >> out
        print >> out, get_rate_matrix_summary(self.max_cheeger_graph)
        print >> out
        print >> out
        print >> out, 'The longest known induced snake-in-the-box path'
        print >> out, 'for a %d-cube has %d states' % (self.ndim, nsnake)
        print >> out, 'with connectivity %f.' % snake_connect
        print >> out
        print >> out, 'Hypercubes have constant algebraic connectivity 2.0'
        print >> out, 'which I hypothesize is an upper bound'
        print >> out, 'of that of their induced subgraphs.'
        return out.getvalue().rstrip()

def get_response_content(fs):
    out = StringIO()
    np.set_printoptions(linewidth=200)
    # look at a bunch of random induced nontrivial connected subgraphs
    f = MyCallback(fs.ndim, fs.pkeep)
    info = combobreaker.run_callable(f, nseconds=5)
    print >> out, info
    return out.getvalue()

