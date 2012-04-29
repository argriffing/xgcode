"""
This is another attempt at what we call neighborhood joining.

I think that neighborhood joining is already a technical term
for something else, but I will use it to mean what I want it to mean.
Here it means a tree reconstruction by divisive clustering.
The divisive clustering is through a supertree object,
whose supervertices are completely connected components
which are graphs of non-super vertices.
The idea is that a single fully connected component is resolved
into a supertree of connected components
each consisting of 2 or 3 vertices.
Both the vertices and the supervertices are labeled by integers.
When a supervertex splits,
both of its child supervertices get new labels.
But the constituent vertices keep their original labels
and a single new constituent vertex is added.
This new constituent vertex becomes a member of both
new supervertices.
Some naming conventions are
v: vertex (an integer index starting from zero)
sv: a supervertex (also an integer index starting from zero)
L: laplacian matrix
"""

import unittest
import itertools
import random
import math

import numpy as np
import scipy
from scipy import linalg

from MatrixUtil import ndot

def harmonic_split_transform(
        sv_big, v_to_svs, sv_to_vs, edge_to_weight,
        v_new, sv_new_a, sv_new_b):
    """
    Transform a supervertex with k>3 vertices to two supervertices.
    The two new supervertices share a single new vertex
    and bipartition the members of the original supervertex.
    This uses the Fiedler vector.
    When the underlying graph is tree-like,
    it splits at exactly the harmonic root of the tree.
    @param sv_big: the supervertex to split
    @param v_to_svs: vertex to set of containing supervertices
    @param sv_to_vs: supervertex to set of contained vertices
    @param edge_to_weight: vertex doubleton to positive edge weight
    @param v_new: new vertex to be
    @param sv_new_a: new supervertex
    @param sv_new_b: new supervertex
    """
    # construct the laplacian matrix
    ord_vs = sorted(sv_to_vs[sv_big])
    n = len(ord_vs)
    L = np.zeros((n, n))
    for i, vi in enumerate(ord_vs):
        for j, vj in enumerate(ord_vs):
            if i != j:
                L[i, j] = -edge_to_weight[frozenset((vi, vj))]
    L -= np.diag(np.sum(L, axis=1))
    # get the fiedler vector
    W, V = scipy.linalg.eigh(L)
    fvect = V.T[1]
    # define the fiedler split of vertices
    a = set(i for i, x in enumerate(fvect) if x < 0)
    b = set(i for i, x in enumerate(fvect) if x >= 0)
    vs_a = sorted(ord_vs[k] for k in a)
    vs_b = sorted(ord_vs[k] for k in b)
    v_to_loading = dict((ord_vs[k], fvect[k]) for k in range(n))
    fa = np.array([v_to_loading[v] for v in vs_a])
    fb = np.array([v_to_loading[v] for v in vs_b])
    # Define a block of the Laplacian matrix
    # chosen according to the signs of the Fiedler vector.
    na = len(vs_a)
    nb = len(vs_b)
    neg_B = np.zeros((na, nb))
    for i, vi in enumerate(vs_a):
        for j, vj in enumerate(vs_b):
            neg_B[i, j] = edge_to_weight[frozenset((vi, vj))]
    # Define the weights to the new vertex using a rank one approximation
    # of a block of the Laplacian matrix.
    U, s, Vh = scipy.linalg.svd(neg_B, full_matrices=False, compute_uv=True)
    u = U.T[0]
    v = Vh[0]
    # SVD is defined only up to the sign of the vectors.
    # Normalize to make them both positive instead of both negative.
    # If they have mixed sign then there is a problem.
    if np.all(u < 0) and np.all(v < 0):
        u *= -1
        v *= -1
    if not (np.all(u > 0) and np.all(v > 0)):
        report_summary(L, neg_B)
        raise ValueError('sign problem with u and v')
    # Decide how to partition the singular value between the two vectors.
    # This is where the harmonicity comes into play.
    sigma = s[0]
    x, y = _get_xy(sigma, u, v, fa, fb)
    z = np.sum(x) + np.sum(y)
    a_block_fail = False
    b_block_fail = False
    for i, j in itertools.combinations(range(na), 2):
        vi = vs_a[i]
        vj = vs_a[j]
        edge = frozenset((vi, vj))
        if edge_to_weight[edge] < x[i] * x[j] / z:
            a_block_fail = True
    for i, j in itertools.combinations(range(nb), 2):
        vi = vs_b[i]
        vj = vs_b[j]
        edge = frozenset((vi, vj))
        if edge_to_weight[edge] < y[i] * y[j] / z:
            b_block_fail = True
    if a_block_fail and b_block_fail:
        report_summary(L, neg_B)
        raise ValueError('inducing unavoidable negative edge weight')
    # Define the set of vertices associated with each new supervertex.
    # Do not add the new vertex yet.
    sv_to_vs[sv_new_a] = set(vs_a)
    sv_to_vs[sv_new_b] = set(vs_b)
    # update the set of supervertices associated with each vertex
    for sv_new in (sv_new_a, sv_new_b):
        for v in sv_to_vs[sv_new]:
            svs = v_to_svs[v]
            v_to_svs[v] = set(sv_new if sv == sv_big else sv for sv in svs)
    # Add the new vertex to both of the new supervertices.
    sv_to_vs[sv_new_a].add(v_new)
    sv_to_vs[sv_new_b].add(v_new)
    v_to_svs[v_new] = set([sv_new_a, sv_new_b])
    # Update the edge weights within each new component.
    for i, j in itertools.combinations(range(na), 2):
        vi = vs_a[i]
        vj = vs_a[j]
        edge = frozenset((vi, vj))
        edge_to_weight[edge] -= x[i] * x[j] / z
    for i, j in itertools.combinations(range(nb), 2):
        vi = vs_b[i]
        vj = vs_b[j]
        edge = frozenset((vi, vj))
        edge_to_weight[edge] -= y[i] * y[j] / z
    # Set the edge weights to the new vertex.
    for xi, vi in zip(x, vs_a):
        edge = frozenset((vi, v_new))
        edge_to_weight[edge] = xi
    for yi, vi in zip(y, vs_b):
        edge = frozenset((vi, v_new))
        edge_to_weight[edge] = yi

def delta_wye_transform(
        sv_big, v_to_svs, sv_to_vs, edge_to_weight,
        v_new, sv_new_a, sv_new_b, sv_new_c):
    """
    Transform a three-vertex supervertex to three two-vertex supervertices.
    """
    # do some error checking
    for v in sv_to_vs[sv_big]:
        if sv_big not in v_to_svs[v]:
            raise ValueError('vertex inclusion error')
    # get the set of vertices
    ord_vs = sorted(sv_to_vs[sv_big])
    # update the weights
    va, vb, vc = ord_vs
    cc = edge_to_weight[frozenset((va, vb))]
    ca = edge_to_weight[frozenset((vb, vc))]
    cb = edge_to_weight[frozenset((vc, va))]
    c1, c2, c3 = _delta_to_wye_conductor(ca, cb, cc)
    edge_to_weight[frozenset((va, v_new))] = c1
    edge_to_weight[frozenset((vb, v_new))] = c2
    edge_to_weight[frozenset((vc, v_new))] = c3
    # update the set of supervertices associated with each vertex
    for v in ord_vs:
        v_to_svs[v].remove(sv_big)
    v_to_svs[va].add(sv_new_a)
    v_to_svs[vb].add(sv_new_b)
    v_to_svs[vc].add(sv_new_c)
    v_to_svs[v_new] = set([sv_new_a, sv_new_b, sv_new_c])
    sv_to_vs[sv_new_a] = set([va, v_new])
    sv_to_vs[sv_new_b] = set([vb, v_new])
    sv_to_vs[sv_new_c] = set([vc, v_new])

def _get_xy(sigma, u, v, f, g):
    """
    This is related to new graph weights when a component is split.
    The x and y are related to edge weights connecting the new vertex
    to existing vertices.
    The returned x and y are scaled u and v.
    @param sigma: singular value
    @param u: singular unit vector
    @param v: singular unit vector
    @param f: fiedler subvector conformant to u
    @param g: fiedler subvector conformant to v
    @return: x, y
    """
    if not np.all(u > 0):
        raise ValueError('input u vector must be positive')
    if not np.all(v > 0):
        raise ValueError('input u vector must be positive')
    if not np.all(f < 0):
        raise ValueError('input f vector must be negative')
    if not np.all(g > 0):
        raise ValueError('input g vector must be positive')
    fu = np.dot(u, f)
    gv = np.dot(v, g)
    sum_u = np.sum(u)
    sum_v = np.sum(v)
    a = sigma * (sum_v - (gv / fu) * sum_u)
    b = sigma * (sum_u - (fu / gv) * sum_v)
    x = a * u
    y = b * v
    return x, y

def _delta_to_wye_resistor(ra, rb, rc):
    """
    Use the wikipedia notation.
    http://upload.wikimedia.org/wikipedia/commons/c/cb/Wye-delta-2.svg
    @return: r1, r2, r3
    """
    r1 = (rb*rc) / (ra + rb + rc)
    r2 = (rc*ra) / (ra + rb + rc)
    r3 = (ra*rb) / (ra + rb + rc)
    return r1, r2, r3

def _delta_to_wye_conductor(ca, cb, cc):
    """
    Analogous to the resistor version.
    Conductor values are reciprocal of resistor values.
    @return: c1, c2, c3
    """
    numerator = cb*cc + cc*ca + ca*cb
    c1 = numerator / ca
    c2 = numerator / cb
    c3 = numerator / cc
    return c1, c2, c3



#FIXME obsolete
def report_summary(L, neg_B):
    U, s, Vh = scipy.linalg.svd(neg_B, full_matrices=False, compute_uv=True)
    x = U.T[0]
    y = Vh[0]
    print 'L:'
    print L
    print
    print 'neg_B:'
    print neg_B
    print
    print 'U:'
    print U
    print
    print 's:'
    print s
    print
    print 'Vh:'
    print Vh
    print 
    print 'U S vh:'
    print ndot(U, np.diag(s), Vh)
    print 
    print 'approx:'
    print s[0] * np.outer(x, y)
    print 
    xycat = x.tolist() + y.tolist()
    signs = set(np.sign(xycat).astype(np.int))
    if set([-1, 1]) <= signs:
        raise ValueError('multiple signs in the concatenated xy')

#FIXME obsolete
def fsplit_random(L):
    nstates = len(L)
    na = random.randrange(1, nstates)
    a = set(random.sample(range(nstates), na))
    b = set(range(nstates)) - a
    return a, b

#FIXME obsolete
def fsplit_fiedler(L):
    nstates = len(L)
    W, V = scipy.linalg.eigh(L)
    fvect = V.T[1]
    a = set(i for i, x in enumerate(fvect) if x < 0)
    b = set(range(nstates)) - a
    return a, b

#FIXME obsolete
def split(
        v_to_svs, sv_to_vs, edge_to_weight, sv_heap,
        v_new, sv_new_a, sv_new_b, fsplit):
    """
    @param v_to_svs: vertex to set of containing supervertices
    @param sv_to_vs: supervertex to set of contained vertices
    @param edge_to_weight: vertex doubleton to positive edge weight
    @param sv_heap: a min heap of (-len(sv_to_vs[sv]), sv) pairs
    @param v_new: new vertex to be
    @param sv_new_a: new supervertex
    @param sv_new_b: new supervertex
    @param fsplit: split the vertices given the laplacian matrix
    """
    # pop the biggest connected component off the heap
    neg_size, sv_big = sv_heap.pop()
    # construct the laplacian matrix
    n = -neg_size
    L = np.zeros((n, n))
    ord_vs = sorted(sv_to_vs[sv_big])
    for i, vi in enumerate(ord_vs):
        for j, vj in enumerate(ord_vs):
            if i != j:
                L[i, j] = -edge_to_weight[frozenset((vi, vj))]
    L -= np.diag(np.sum(L, axis=1))
    # define the split using the supplied split function
    a, b = fsplit(L)
    vs_a = sorted(ord_vs[k] for k in a)
    vs_b = sorted(ord_vs[k] for k in b)
    # Define a block of the Laplacian matrix
    # chosen according to the signs of the Fiedler vector.
    na = len(vs_a)
    nb = len(vs_b)
    neg_B = np.zeros((na, nb))
    for i, vi in enumerate(vs_a):
        for j, vj in enumerate(vs_b):
            neg_B[i, j] = edge_to_weight[frozenset((vi, vj))]
    # Define the weights to the new vertex using a rank one approximation
    # of a block of the Laplacian matrix.
    U, s, Vh = scipy.linalg.svd(neg_B, full_matrices=False, compute_uv=True)
    X = U.T[0] * math.sqrt(s[0])
    Y = Vh[0] * math.sqrt(s[0])
    if np.all(X < 0) and np.all(Y < 0):
        X *= -1
        Y *= -1
    if not (np.all(X > 0) and np.all(Y > 0)):
        report_summary(L, neg_B)
        raise ValueError('sign problem with X and Y')
    a_block_fail = False
    b_block_fail = False
    for i, j in itertools.combinations(range(na), 2):
        vi = vs_a[i]
        vj = vs_a[j]
        edge = frozenset((vi, vj))
        if edge_to_weight[edge] < X[i] * X[j]:
            a_block_fail = True
    for i, j in itertools.combinations(range(nb), 2):
        vi = vs_b[i]
        vj = vs_b[j]
        edge = frozenset((vi, vj))
        if edge_to_weight[edge] < Y[i] * Y[j]:
            b_block_fail = True
    if a_block_fail and b_block_fail:
        report_summary(L, neg_B)
        raise ValueError('inducing unavoidable negative edge weight')
    """
    #approx = s[0] * np.outer(x, y)
    #sign_observed = np.sign(approx).astype(np.int)
    #sign_expected = -1 * np.ones((na, nb))
    #if not np.array_equal(sign_observed, sign_expected):
        #raise ValueError('rank one approximation is not negative')
    # define the set of vertices associated with each new supervertex
    sv_to_vs[sv_new_a] = set(v for v, x in zip(ord_vs, fvect) if x < 0)
    sv_to_vs[sv_new_b] = set(v for v, x in zip(ord_vs, fvect) if x >= 0)
    # update the v to svs dict
    for sv_new in (sv_new_a, sv_new_b):
        for v in sv_to_vs[sv_new]:
            svs = v_to_svs[v]
            v_to_svs[v] = set(sv_new_a if sv == sv_big else sv for sv in svs)
    """

#FIXME obsolete
def harmonic_split(
        v_to_svs, sv_to_vs, edge_to_weight, sv_heap,
        v_new, sv_new_a, sv_new_b):
    """
    This uses the Fiedler vector.
    When the underlying graph is tree-like,
    it splits at exactly the harmonic root of the tree.
    @param v_to_svs: vertex to set of containing supervertices
    @param sv_to_vs: supervertex to set of contained vertices
    @param edge_to_weight: vertex doubleton to positive edge weight
    @param sv_heap: a min heap of (-len(sv_to_vs[sv]), sv) pairs
    @param v_new: new vertex to be
    @param sv_new_a: new supervertex
    @param sv_new_b: new supervertex
    """
    # pop the biggest connected component off the heap
    neg_size, sv_big = sv_heap.pop()
    # construct the laplacian matrix
    n = -neg_size
    L = np.zeros((n, n))
    ord_vs = sorted(sv_to_vs[sv_big])
    for i, vi in enumerate(ord_vs):
        for j, vj in enumerate(ord_vs):
            if i != j:
                L[i, j] = -edge_to_weight[frozenset((vi, vj))]
    L -= np.diag(np.sum(L, axis=1))
    # define Fiedler split
    W, V = scipy.linalg.eigh(L)
    fvect = V.T[1]
    a = set(i for i, x in enumerate(fvect) if x < 0)
    b = set(range(nstates)) - a
    vs_a = sorted(ord_vs[k] for k in a)
    vs_b = sorted(ord_vs[k] for k in b)
    v_to_loading = dict((ord_vs[k], fvect[k]) for k in range(n))
    fa = np.array([v_to_loading[v] for v in vs_a])
    fb = np.array([v_to_loading[v] for v in vs_b])
    # Define a block of the Laplacian matrix
    # chosen according to the signs of the Fiedler vector.
    na = len(vs_a)
    nb = len(vs_b)
    neg_B = np.zeros((na, nb))
    for i, vi in enumerate(vs_a):
        for j, vj in enumerate(vs_b):
            neg_B[i, j] = edge_to_weight[frozenset((vi, vj))]
    # Define the weights to the new vertex using a rank one approximation
    # of a block of the Laplacian matrix.
    U, s, Vh = scipy.linalg.svd(neg_B, full_matrices=False, compute_uv=True)
    u = U.T[0]
    v = Vh[0]
    # SVD is defined only up to the sign of the vectors.
    # Normalize to make them both positive instead of both negative.
    # If they have mixed sign then there is a problem.
    if np.all(u < 0) and np.all(v < 0):
        u *= -1
        v *= -1
    if not (np.all(u > 0) and np.all(v > 0)):
        report_summary(L, neg_B)
        raise ValueError('sign problem with u and v')
    # Decide how to partition the singular value between the two vectors.
    # This is where the harmonicity comes into play.
    alpha = math.sqrt(-s * np.dot(u, fa) / np.dot(v, fb))
    beta = s / alpha
    x = u * alpha
    y = u * beta
    a_block_fail = False
    b_block_fail = False
    for i, j in itertools.combinations(range(na), 2):
        vi = vs_a[i]
        vj = vs_a[j]
        edge = frozenset((vi, vj))
        if edge_to_weight[edge] < x[i] * x[j]:
            a_block_fail = True
    for i, j in itertools.combinations(range(nb), 2):
        vi = vs_b[i]
        vj = vs_b[j]
        edge = frozenset((vi, vj))
        if edge_to_weight[edge] < y[i] * y[j]:
            b_block_fail = True
    if a_block_fail and b_block_fail:
        report_summary(L, neg_B)
        raise ValueError('inducing unavoidable negative edge weight')


class TestNHJ(unittest.TestCase):

    def test_delta_to_wye_resistor(self):
        """
        Use an example from the internet.
        http://www.tina.com/English/tina/course/5wye/wye.htm
        """
        ra, rb, rc = 14.0, 7.0, 3.5
        expected_triple = (1.0, 2.0, 4.0)
        observed_triple = delta_to_wye_resistor(ra, rb, rc)
        self.assertTrue(np.allclose(expected_triple, observed_triple))

    def test_delta_to_wye_conductor(self):
        ca, cb, cc = 1 / 14.0, 1 / 7.0, 1 / 3.5
        expected_triple = (1 / 1.0, 1 / 2.0, 1 / 4.0)
        observed_triple = delta_to_wye_conductor(ca, cb, cc)
        self.assertTrue(np.allclose(expected_triple, observed_triple))



if __name__ == '__main__':
    unittest.main()

