"""
Visualize steps of an exact spectral reconstruction of a tree.
"""

from StringIO import StringIO
import itertools
import math

import numpy as np
import scipy
from scipy import linalg
import pydot

import Form
import FormOut
import nhj
import Util
from MatrixUtil import ndot

g_upper_adjacency = np.array([
    [0, 0, 0, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 0, 2, 0, 0],
    [0, 0, 0, 0, 0, 0, 3, 0],
    [0, 0, 0, 0, 0, 0, 2, 0],
    [0, 0, 0, 0, 0, 0, 0, 1],
    [0, 0, 0, 0, 0, 0, 0, 1],
    [0, 0, 0, 0, 0, 0, 0, 1],
    [0, 0, 0, 0, 0, 0, 0, 0]])

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            ]
    return form_objects

def get_form_out():
    return FormOut.Html('reconstruction')

"""
def foo():
    n = self.ntaxa
    # sample a laplacian matrix
    L = sample_laplacian_matrix(n)
    # get the new vertex and the two new replacement supervertices
    v_gen = itertools.count(n)
    sv_gen = itertools.count(1)
    v_new = next(v_gen)
    sv_new_a = next(sv_gen)
    sv_new_b = next(sv_gen)
    # define the initial structures
    v_to_svs = dict((v, set([0])) for v in range(n))
    sv_to_vs = {0 : set(range(n))}
    edge_to_weight = {}
    for pair in itertools.combinations(range(n), 2):
        edge_to_weight[frozenset(pair)] = random.expovariate(1)
    sv_heap = [(-n, 0)]
    # do the split
    nhj.split(
            v_to_svs, sv_to_vs, edge_to_weight, sv_heap,
            v_new, sv_new_a, sv_new_b, self.fsplit)
"""

def get_svg(active_svs, sv_to_vs, v_to_name, v_to_svs, edge_to_weight):
    # get the set of active vertices
    vs = set()
    for sv in active_svs:
        vs.update(sv_to_vs[sv])
    # for debugging show the vertices and the supervertices
    #print vs
    # initialize the graph
    pydot_graph = pydot.Dot(
            graph_type='graph',
            overlap='0',
            sep='0.01',
            size='8, 8',
            )
    """
    pydot_graph = pydot.Dot(
            graph_type='digraph',
            size='%s, 8' % width_inches,
            rankdir=rankdir,
            ratio='compress')
    """
    # define the pydot node objects with the right names
    v_to_pydot_node = dict((v, pydot.Node(v_to_name[v])) for v in vs)
    # define the edges
    edge_to_pydot_edge = {}
    for sv in active_svs:
        nvertices = len(sv_to_vs[sv])
        for pair in itertools.combinations(sv_to_vs[sv], 2):
            edge = frozenset(pair)
            distance = 1 / edge_to_weight[edge]
            va, vb = pair
            pna = v_to_pydot_node[va]
            pnb = v_to_pydot_node[vb]
            if nvertices == 2:
                # annotate the edge with the branch length
                label = '%.3f' % distance
                pydot_edge = pydot.Edge(pna, pnb, label=label)
            else:
                # do not annotate the edge with the branch length
                pydot_edge = pydot.Edge(pna, pnb)
            edge_to_pydot_edge[edge] = pydot_edge
    # add the nodes and edges
    for pydot_node in v_to_pydot_node.values():
        pydot_graph.add_node(pydot_node)
    for pydot_edge in edge_to_pydot_edge.values():
        pydot_graph.add_edge(pydot_edge)
    # do the physical layout and create the svg string
    tmp_path = Util.create_tmp_file(data=None, prefix='tmp', suffix='.svg')
    pydot_graph.write_svg(tmp_path, prog='neato')
    with open(tmp_path) as fin:
        svg_str = fin.read()
    # return the svg except for the first few lines
    svg_str = '\n'.join(svg_str.splitlines()[6:])
    return svg_str

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
    rc = 1 / edge_to_weight[frozenset((va, vb))]
    ra = 1 / edge_to_weight[frozenset((vb, vc))]
    rb = 1 / edge_to_weight[frozenset((vc, va))]
    r1 = (rb*rc) / (ra + rb + rc)
    r2 = (rc*ra) / (ra + rb + rc)
    r3 = (ra*rb) / (ra + rb + rc)
    edge_to_weight[frozenset((va, v_new))] = 1 / r1
    edge_to_weight[frozenset((vb, v_new))] = 1 / r2
    edge_to_weight[frozenset((vc, v_new))] = 1 / r3
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

def get_xy(sigma, u, v, f, g):
    """
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
    x, y = get_xy(sigma, u, v, fa, fb)
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


def get_response_content(fs):
    A = g_upper_adjacency + g_upper_adjacency.T
    L = (np.diag(np.sum(A, axis=1)) - A).astype(np.float64)
    # define the number of pendant vertices (leaves)
    p = 5
    # define the fully connected schur complement graph as a Laplacian matrix
    G = L[:p,:p] - ndot(L[:p,p:], np.linalg.pinv(L[p:,p:]), L[p:,:p])
    # init the tree reconstruction state
    v_to_name = dict((v, 'P' + chr(ord('a') + v)) for v in range(p))
    v_to_svs = dict((v, set([0])) for v in range(p))
    sv_to_vs = {0 : set(range(p))}
    edge_to_weight = {}
    for pair in itertools.combinations(range(p), 2):
        edge_to_weight[frozenset(pair)] = -G[pair]
    # pairs like (-(number of vertices in supervertex sv), supervertex sv)
    active_svs = set([0])
    # initialize the sources of unique vertex and supervertex identifiers
    v_gen = itertools.count(5)
    sv_gen = itertools.count(1)
    # write the output
    out = StringIO()
    print >> out, '<html>'
    print >> out, '<body>'
    for count_pos in itertools.count(1):
        # add the graph rendering before the decomposition at this stage
        print >> out, '<div>'
        print >> out, get_svg(
                active_svs, sv_to_vs, v_to_name, v_to_svs, edge_to_weight)
        print >> out, '</div>'
        # update the splits
        next_active_svs = set()
        # svs can be decomposed independently in arbitrary order
        alpha_index_gen = itertools.count()
        for sv in active_svs:
            if len(sv_to_vs[sv]) > 2:
                v_new = next(v_gen)
                sv_new_a = next(sv_gen)
                sv_new_b = next(sv_gen)
                alpha_index = next(alpha_index_gen)
                alpha = chr(ord('a') + alpha_index)
                v_to_name[v_new] = 'R%s%s' % (count_pos, alpha)
                next_active_svs.add(sv_new_a)
                next_active_svs.add(sv_new_b)
                if len(sv_to_vs[sv]) == 3:
                    sv_new_c = next(sv_gen)
                    delta_wye_transform(
                            sv, v_to_svs, sv_to_vs, edge_to_weight,
                            v_new, sv_new_a, sv_new_b, sv_new_c)
                    next_active_svs.add(sv_new_c)
                else:
                    harmonic_split_transform(
                            sv, v_to_svs, sv_to_vs, edge_to_weight,
                            v_new, sv_new_a, sv_new_b)
            else:
                next_active_svs.add(sv)
        # if the set of active svs has not changed then we are done
        if active_svs == next_active_svs:
            break
        else:
            active_svs = next_active_svs
    print >> out, '</html>'
    print >> out, '</body>'
    return out.getvalue()

