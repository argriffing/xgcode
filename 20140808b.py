"""
Check properties of Fiedler-analog vector of normalized Laplacian.

In particular check if the center-of-mass cut of the Fiedler-analog vector
of the normalized Laplacian matrix of the Schur complement tree graph
always minimizes the sum of within-cluster variance of the two clusters
induced by the cut.

"""
from StringIO import StringIO
import random
import math

import networkx as nx
import numpy as np
import scipy.linalg
from scipy.linalg import inv

import Form
import FormOut
import const
import Euclid
import combobreaker


def get_form():
    """
    @return: a list of form objects
    """
    # define the form objects
    form_objects = [
            ]
    return form_objects


def get_form_out():
    return FormOut.Report()


def get_best_sum_of_variance_split(v):
    # input should be a sorted vector
    # return (best sum of variances, nfirst)
    n = len(v)
    pairs = []
    for nfirst in range(1, n):
        value = np.var(v[:nfirst]) + np.var(v[nfirst:])
        pairs.append((value, nfirst))
    return min(pairs)


def get_center_of_mass_split(v):
    # input should be a sorted vector
    # return nfirst
    n = len(v)
    for nfirst in range(1, n):
        if v[nfirst-1] * v[nfirst] < 0:
            return nfirst
    raise Exception


class State(object):
    def __init__(self):
        # sample a random rooted tree with weighted edges
        nedges_min, nedges_max = 3, 20
        nedges = random.randrange(nedges_min, nedges_max+1)
        T = nx.Graph()
        root = 0
        T.add_node(root)
        for edge_idx in range(nedges):
            nb = edge_idx+1
            na = random.randrange(nb)
            weight = random.expovariate(1)
            T.add_edge(na, nb, weight=weight)

        # define the node list including leaves and internal vertices
        degs = T.degree()
        leaves = [v for v in T if degs[v] == 1]
        internals = [v for v in T if degs[v] > 1]
        nodes = leaves + internals

        # get the laplacian matrix
        L = nx.linalg.laplacianmatrix.laplacian_matrix(
                T, nodes, weight='weight').A
        if len(L.shape) != 2 or L.shape[0] != L.shape[1]:
            raise Exception(str(L.shape))

        # compute the schur complement of the laplacian matrix
        node_partition_string = str((leaves, internals))
        if len(leaves) < 2 or len(internals) < 2:
            raise Exception(node_partition_string)
        n = len(leaves)
        try:
            L_schur = L[:n, :n] - L[:n, n:].dot(inv(L[n:, n:])).dot(L[n:, :n])
        except ValueError as e:
            raise Exception((node_partition_string, L[n:, n:]))

        # get the normalized laplacian
        m = np.reciprocal(np.sqrt(np.diag(L_schur)))
        L_norm = L_schur * np.outer(m, m)

        # get the fiedler vector from the schur complement matrix
        w, V = scipy.linalg.eigh(L_norm)
        v = V[:, 1]

        # store some stuff
        self.T = T
        self.leaves = leaves
        self.internals = internals
        self.nodes = nodes
        self.L = L
        self.L_schur = L_schur
        self.L_norm = L_norm
        self.w = w
        self.V = V
        self.v = v

    def is_counterexample(self):
        # check each split for the min variance and the center of mass
        self.vsorted = sorted(self.v)
        info = get_best_sum_of_variance_split(self.vsorted)
        self.sum_of_variances, self.var_nfirst = info
        self.fiedler_nfirst = get_center_of_mass_split(self.vsorted)
        if self.var_nfirst != self.fiedler_nfirst:
            return True

    def __str__(self):
        out = StringIO()
        print >> out, 'weights of edges of tree:'
        for edge in self.T.edges():
            va, vb = edge
            weight = self.T[va][vb]['weight']
            print >> out, va, vb, weight
        print >> out
        print >> out, 'full laplacian matrix:'
        print >> out, self.L
        print >> out
        print >> out, 'schur complement laplacian matrix:'
        print >> out, self.L_schur
        print >> out
        print >> out, 'normalized schur complement laplacian matrix:'
        print >> out, self.L_norm
        print >> out
        print >> out, 'fiedler vector:'
        print >> out, self.v
        print >> out
        print >> out, 'sorted fiedler vector:'
        print >> out, self.vsorted
        print >> out
        print >> out, 'min sum of-variance size of |A|:', self.var_nfirst
        print >> out, 'center of mass size of |A|:', self.fiedler_nfirst
        return out.getvalue().rstrip()

def get_response_content(fs):

    # set up print options
    np.set_printoptions(
            linewidth=1000000,
            threshold=1000000,
            )
    out = StringIO()

    ntries = 100
    found = False
    for i in range(ntries):
        state = State()
        if state.is_counterexample():
            found = True
            print >> out, 'found a counterexample on attempt', i+1
            print >> out
            print >> out, state
            break

    if not found:
        print >> out, 'gave up after', ntries, 'tries'

    # show the result
    return out.getvalue()

