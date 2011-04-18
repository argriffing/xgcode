"""
This is a function oriented rather than object oriented approach to trees.

The basic trees in this module are unrooted but have labeled vertices.
A vertex is something hashable.
Vertices of a tree are distinct.
Edges are frozensets of size two with vertex elements.
The tree topology is defined by a set of edges.
Branch lengths are defined by a dictionary
where dictionary keys are edges
and dictionary values are floating point branch lengths.
A directed topology is a set of vertex pairs.
The order of vertices is sometimes used to facilitate testing.
Sometimes the following notations will be used: T is a set of undirected edges,
B is a map from vertex doubletons to branch lengths,
R is a set of directed edges, D is a tree distance matrix,
L is a tree Laplacian matrix.
Matrices are numpy arrays.
"""

import unittest
from StringIO import StringIO
from collections import defaultdict

import numpy as np

import MatrixUtil


def mkedge(a, b):
    """
    @param a: a vertex
    @param b: another vertex
    @return: the unordered edge between the two vertices
    """
    return frozenset([a, b])


# In this section we define functions on directed trees.

def get_v_to_sinks(R):
    """
    @param R: a directed topology
    @return: a map from a vertex to a set of sinks
    """
    d = defaultdict(set)
    for a, b in R:
        d[a].add(b)
    return d

def get_v_to_source(R):
    """
    @param R: a directed topology
    @return: a map from a vertex to its source
    """
    return dict((b, a) for a, b in R)

def get_R_vertices(R):
    """
    @param R: a directed topology
    @return: the set of vertices
    """
    V = set()
    for a, b in R:
        V.add(a)
        V.add(b)
    return V

def get_root(R):
    """
    @param R: a directed topology
    @return: the root vertex
    """
    root_set = get_R_vertices(R) - set(get_v_to_source(R))
    return sorted(root_set)[0]

def _v_to_newick(v_to_sinks, v):
    """
    This is part of writing a pure topology newick string.
    Note that this function does not add the semicolon termination.
    @param v_to_sinks: a map from a vertex to a set of sinks
    @param v: a subtree root vertex
    @return: a chunk of a newick string
    """
    if v not in v_to_sinks:
        return str(v)
    sinks = sorted(v_to_sinks[v])
    arr = [_v_to_newick(v_to_sinks, x) for x in sinks]
    return '(' + ', '.join(arr) + ')' + str(v)

def R_to_newick(R):
    """
    @param R: a directed topology
    @return: a newick string
    """
    r = get_root(R)
    return _v_to_newick(get_v_to_sinks(R), r) + ';'

def _Bv_to_newick(v_to_source, v_to_sinks, B, v):
    """
    This is part of writing a newick string with branch lengths.
    Note that this function does not add the semicolon termination.
    @param v_to_source: a map from a vertex to its source
    @param v_to_sinks: a map from a vertex to a set of sinks
    @param B: branch lengths
    @param v: a subtree root vertex
    @return: a chunk of a newick string
    """
    if v in v_to_source:
        # a vertex that has a source should record its distance to its source
        blen = B[mkedge(v, v_to_source[v])]
        suffix = ':' + str(blen)
    else:
        suffix = ''
    if v not in v_to_sinks:
        return str(v) + suffix
    sinks = sorted(v_to_sinks[v])
    arr = [_Bv_to_newick(v_to_source, v_to_sinks, B, x) for x in sinks]
    return '(' + ', '.join(arr) + ')' + str(v) + suffix

def RB_to_newick(R, B):
    """
    @param R: a directed topology
    @param B: branch lengths
    @return: a newick string
    """
    r = get_root(R)
    return _Bv_to_newick(get_v_to_source(R), get_v_to_sinks(R), B, r) + ';'

def R_to_preorder(R):
    """
    Return the preorder list of vertices.
    This list has the property that every child vertex appears later
    in the list than its parent vertex.
    @param R: a directed topology
    @return: the preorder list of vertices
    """
    v_to_sinks = get_v_to_sinks(R)
    r = get_root(R)
    shell = set([r])
    visited = set()
    arr = []
    while shell:
        next_shell = set()
        for v in sorted(shell):
            arr.append(v)
            for adj in sorted(v_to_sinks[v]):
                next_shell.add(adj)
        shell = next_shell
    return arr

def RB_to_D(R, B, vertices):
    """
    Get the matrix of distances among the given vertices.
    @param R: a directed topology
    @param B: branch lengths
    @param vertices: ordered vertices
    @return: distance matrix
    """
    v_to_source = get_v_to_source(R)
    # An intermediate representation is
    # a dense map from an edge to a distance.
    visited = set()
    edge_to_distance = {}
    for v in R_to_preorder(R):
        if v in v_to_source:
            parent = v_to_source[v]
            excess = B[mkedge(v, parent)]
            edge_to_distance[mkedge(v, parent)] = excess
            for n in visited:
                if n != parent:
                    distance = edge_to_distance[mkedge(parent, n)] + excess
                    edge_to_distance[mkedge(v,n)] = distance
        visited.add(v)
    # Convert this intermediate representation
    # to the requested distance matrix.
    N = len(vertices)
    v_to_index = dict((v, i) for i, v in enumerate(vertices))
    D = np.zeros((N, N))
    for (a, b), d in edge_to_distance.items():
        if (a in v_to_index) and (b in v_to_index):
            i = v_to_index[a]
            j = v_to_index[b]
            D[i,j] = d
            D[j,i] = d
    return D


# In this section we define functions on undirected trees.

def T_to_order(T):
    """
    Get an ordered sequence of vertices.
    The order should have the property that every graph induced by
    an initial subsequence should be connected.
    @return: a list of ordered vertices
    """
    return R_to_preorder(T_to_R_canonical(T))

def TB_to_D(T, B, vertices):
    """
    Get the matrix of distances among the given vertices.
    @param T: topology
    @param B: branch lengths
    @param vertices: ordered vertices
    @return: distance matrix
    """
    return RB_to_D(T_to_R_canonical(T), B, vertices)

def TB_to_L_block(T, B, row_vertices, col_vertices):
    """
    Get the Laplacian matrix.
    @param T: topology
    @param B: branch lengths
    @param row_vertices: ordered vertices corresponding to rows
    @param col_vertices: ordered vertices corresponding to columns
    @return: a block of a laplacian matrix
    """
    L_part = np.zeros((len(row_vertices), len(col_vertices)))
    v_to_degree = TB_to_v_to_degree(T, B)
    for i, v_row in row_vertices:
        for j, v_col in col_vertices:
            if v_row == v_col:
                L_part[i, j] = v_to_degree[v_row]
            else:
                edge = mkedge(v_row, v_col)
                if edge in B:
                    w = 1.0 / B[edge]
                    L_part[i, j] = -w
    return L_part

def TB_to_L_principal(T, B, vertices):
    """
    Get a principal submatrix of the Laplacian matrix.
    @param T: topology
    @param B: branch lengths
    @param vertices: ordered vertices
    @return: a block of a laplacian matrix
    """
    return TB_to_L_block(T, B, vertices, vertices)

def T_to_v_to_neighbors(T):
    """
    Get the map from a vertex to the set of its adjacent vertices.
    @param T: topology
    @return: a map from a vertex to the set of its adjacent vertices.
    """
    d = defaultdict(set)
    for a, b in T:
        d[a].add(b)
        d[b].add(a)
    return d

def TB_to_v_to_degree(T, B):
    """
    Get the map from a vertex to the undirected weighted vertex degree.
    @param T: topology
    @param B: branch lengths
    @return: a map from a vertex to sum of weights to adjacent vertices
    """
    d = defaultdict(float)
    for edge in T:
        w = 1.0 / B[edge]
        a, b = edge
        d[a] += w
        d[b] += w
    return d

def T_to_v_to_degree(T):
    """
    Get the map from a vertex to the undirected unweighted vertex degree.
    @param T: topology
    @return: a map from a vertex to the size of its set of adjacent vertices
    """
    d = defaultdict(int)
    for a, b in T:
        d[a] += 1
        d[b] += 1
    return d

def T_to_root(T):
    """
    Pick an arbitrary root vertex given an undirected topology.
    No other vertex will be adjacent to more vertices.
    Of the vertices adjacent to the same number of vertices,
    no other vertex will precede it in the natural ordering of the vertices.
    """
    v_to_degree = T_to_v_to_degree(T)
    pairs = [(-deg, v) for v, deg in v_to_degree.items()]
    best_neg_deg, best_v = min(pairs)
    return best_v

def T_to_leaves(T):
    """
    @return: sorted degree one vertices
    """
    v_to_degree = T_to_v_to_degree(T)
    return sorted(v for v, degree in v_to_degree.items() if degree == 1)

def T_to_internal_vertices(T):
    """
    @return: sorted vertices of degree greater than one
    """
    v_to_degree = T_to_v_to_degree(T)
    return sorted(v for v, degree in v_to_degree.items() if degree > 1)

def T_to_R_specific(T, r):
    """
    Convert an unrooted topology to a directed topology.
    Use a root specified by the caller.
    @param T: topology
    @return: a directed topology
    """
    R = set()
    v_to_neighbors = T_to_v_to_neighbors(T)
    visited = set()
    shell = set([r])
    while shell:
        next_shell = set()
        for v in shell:
            for adj in v_to_neighbors[v]:
                if adj not in visited:
                    R.add((v, adj))
                    next_shell.add(adj)
        visited |= shell
        shell = next_shell
    return R

def T_to_R_canonical(T):
    """
    Convert an unrooted topology to a directed topology.
    Use a canonical root.
    @param T: topology
    @return: a directed topology
    """
    return T_to_R_specific(T, T_to_root(T))

def T_to_newick(T):
    """
    Get a newick string from an unweighted topology.
    @param T: topology
    @return: newick string
    """
    return R_to_newick(T_to_R_canonical(T))

def TB_to_newick(T, B):
    """
    Get a newick string from a weighted topology.
    @param T: topology
    @param B: branch lengths
    @return: newick string
    """
    return RB_to_newick(T_to_R_canonical(T), B)

g_example_T = set([
    mkedge(2,1),
    mkedge(2,3),
    mkedge(2,4),
    mkedge(3,5),
    mkedge(3,6),
    mkedge(6,7)])

g_example_B = {
        mkedge(2,1) : 1,
        mkedge(2,3) : 2,
        mkedge(2,4) : 2,
        mkedge(3,5) : 3,
        mkedge(3,6) : 3,
        mkedge(6,7) : 3}

class TestFtree(unittest.TestCase):

    def test_leaves(self):
        observed = T_to_leaves(g_example_T)
        expected = [1,4,5,7]
        self.assertEqual(observed, expected)

    def test_internal_vertices(self):
        observed = T_to_internal_vertices(g_example_T)
        expected = [2,3,6]
        self.assertEqual(observed, expected)

    def test_unweighted_newick(self):
        observed = T_to_newick(g_example_T)
        expected = '(1, (5, (7)6)3, 4)2;'
        self.assertEqual(observed, expected)

    def test_weighted_newick(self):
        observed = TB_to_newick(g_example_T, g_example_B)
        expected = '(1:1, (5:3, (7:3)6:3)3:2, 4:2)2;'
        self.assertEqual(observed, expected)

    def test_distance_matrix(self):
        vertices = (3,4,5)
        observed = TB_to_D(g_example_T, g_example_B, vertices)
        expected = np.array([
            [0, 4, 3],
            [4, 0, 7],
            [3, 7, 0]], dtype=float)
        self.assertTrue(np.allclose(observed, expected))


if __name__ == '__main__':
    unittest.main()

