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
The functions that read the Newick strings
do not save enough information to reconstruct the input strings.
"""

import unittest
from StringIO import StringIO
from collections import defaultdict

import numpy as np
import scipy.linalg

import MatrixUtil


# In this section we define helper functions.

def invseq(seq):
    """
    This is a helper function.
    @param seq: an iterable
    @return: a map from the value to the sequence index
    """
    return dict((v, i) for i, v in enumerate(seq))

def mkedge(a, b):
    """
    This is a helper function.
    @param a: a vertex
    @param b: another vertex
    @return: the unordered edge between the two vertices
    """
    return frozenset([a, b])


# In this section we define validation functions.

def TB_assert_branch_lengths(T, B):
    extra = set(B) - set(T)
    missing = set(T) - set(B)
    if not B:
        msg = 'no branch lengths were found'
        raise ValueError(msg)
    if extra:
        msg = 'lengths are specified for nonexistent branches'
        raise ValueError(msg)
    if missing:
        msg = 'some branches do not have lengths'
        raise ValueError(msg)

def RB_assert_branch_lengths(R, B):
    TB_assert_branch_lengths(R_to_T(B), B)


# In this section we define functions on directed trees.

def R_to_v_to_sinks(R):
    """
    @param R: a directed topology
    @return: a map from a vertex to a set of sinks
    """
    d = defaultdict(set)
    for a, b in R:
        d[a].add(b)
    return d

def R_to_v_to_source(R):
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

def R_to_root(R):
    """
    @param R: a directed topology
    @return: the root vertex
    """
    root_set = get_R_vertices(R) - set(R_to_v_to_source(R))
    return sorted(root_set)[0]

def R_to_T(R):
    """
    @param R: a directed topology
    @return: an undirected topology
    """
    return set(frozenset(d_edge) for d_edge in R)

def R_to_postorder(R):
    """
    Return the postorder list of vertices.
    This list has the property that every child vertex appears earlier
    in the list than its parent vertex.
    @param R: a directed topology
    @return: the postorder list of vertices
    """
    return list(reversed(R_to_preorder(R)))

def R_to_preorder(R):
    """
    Return the preorder list of vertices.
    This list has the property that every child vertex appears later
    in the list than its parent vertex.
    @param R: a directed topology
    @return: the preorder list of vertices
    """
    return R_to_subtree_preorder(R, R_to_root(R))

def R_to_vertex_partition(R):
    """
    @param R: a directed topology
    @return: set of frozensets of vertices
    """
    v_to_sinks = R_to_v_to_sinks(R)
    roots = v_to_sinks[R_to_root(R)]
    return set(frozenset(R_to_subtree_preorder(R, r)) for r in roots)

def R_to_subtree_preorder(R, r):
    """
    @param R: directed topology
    @param r: subtree root
    @return: list of vertices
    """
    v_to_sinks = R_to_v_to_sinks(R)
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
    v_to_source = R_to_v_to_source(R)
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
    v_to_index = invseq(vertices)
    D = np.zeros((N, N))
    for (a, b), d in edge_to_distance.items():
        if (a in v_to_index) and (b in v_to_index):
            i = v_to_index[a]
            j = v_to_index[b]
            D[i,j] = d
            D[j,i] = d
    return D


# In this section we define functions on undirected trees.

def T_to_v_to_centrality(T):
    """
    Get a certain kind of centrality for each vertex.
    The particular kind of centerality for this function
    is the number of times that you have to
    strip off the terminal branches of the tree before
    the vertex of interest has degree one.
    So a leaf of the tree has zero centrality.
    """
    v_to_cent = {}
    v_to_neighbors = T_to_v_to_neighbors(T)
    cent = 0
    # The shell is the set of leaves at the current centrality.
    shell = set(T_to_leaves(T))
    while shell:
        # mark the centrality for the given shell
        for v in shell:
            v_to_cent[v] = cent
        cent += 1
        # Get all vertices which are candidates for the next shell.
        next_possible_shell = set()
        for v in shell:
            for adj in v_to_neighbors[v]:
                if adj not in v_to_cent:
                    next_possible_shell.add(adj)
        # Find the next shell vertices adjacent to at most
        # a single vertex whose centrality is unknown.
        shell = set()
        for p in next_possible_shell:
            nfree = sum(1 for a in v_to_neighbors[p] if a not in v_to_cent)
            if nfree < 2:
                shell.add(p)
    return v_to_cent

def T_to_edges(T):
    """
    @param T: a topology
    @return: a list of undirected edges
    """
    return [frozenset(d_edge) for d_edge in T_to_outside_in_edges(T)]

def T_to_outside_in_edges(T):
    """
    @param T: a topology
    @return: an ordered list of directed edges
    """
    edges = []
    vertices = T_to_outside_in(T)
    v_to_cent = T_to_v_to_centrality(T)
    v_to_neighbors = T_to_v_to_neighbors(T)
    for v in vertices[:-1]:
        max_cent, best_adj = max((v_to_cent[x], x) for x in v_to_neighbors[v])
        edges.append((v, best_adj))
    return edges

def T_to_incidence_matrix(T, vertices):
    """
    The order and signs of edges is defined by the outside in edges function.
    @param T: a topology
    @param vertices: ordered vertices
    @return: an unweighted incidence matrix
    """
    edges = T_to_outside_in_edges(T)
    nrows = len(edges) + 1
    ncols = len(edges)
    v_to_index = invseq(vertices)
    C = np.zeros((nrows, ncols))
    for i, (x, y) in enumerate(edges):
        a = v_to_index[x]
        b = v_to_index[y]
        C[a, i] = 1
        C[b, i] = -1
    return C

def TB_to_weight_matrix(T, B):
    """
    The order of edges is defined by the outside in edges function.
    @param T: a topology
    @param B: branch lengths
    @return: a diagonal weight matrix
    """
    edges = T_to_outside_in_edges(T)
    w = [1.0 / B[frozenset(e)] for e in edges]
    return np.diag(w)

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
    for i, v_row in enumerate(row_vertices):
        for j, v_col in enumerate(col_vertices):
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

def TB_to_G(T, B, vertices):
    """
    Get the Gower matrix.
    Note that this should be the pseudoinverse of L_schur.
    @param T: topology
    @param B: branch lengths
    @param vertices: ordered vertices
    @return: the Gower matrix
    """
    D = TB_to_D(T, B, vertices)
    HDH = MatrixUtil.double_centered(D)
    return -0.5 * HDH

def TB_to_L_schur(T, B, vertices):
    """
    Get a Schur complement in the Laplacian matrix.
    The specified vertices are the ones to keep.
    So the returned matrix will by a Schur complement Laplacian matrix
    which is square symmetric singular positive semidefinite
    and whose number of rows and columns
    is the same as the number of provided vertices.
    The Schur complement is in the full Laplacian
    defined by all vertices in edges of T.
    Note that this should be the pseudoinverse of G.
    @param T: topology
    @param B: branch lengths
    @param vertices: ordered vertices
    @return: a Schur complement matrix
    """
    # define the list of removed vertices
    removed  = sorted(set(T_to_order(T)) - set(vertices))
    Laa = TB_to_L_block(T, B, vertices, vertices)
    Lab = TB_to_L_block(T, B, vertices, removed)
    Lba = TB_to_L_block(T, B, removed, vertices)
    Lbb = TB_to_L_block(T, B, removed, removed)
    # use the Schur complement definition directly
    L_schur = Laa - np.dot(Lab, np.dot(np.linalg.pinv(Lbb), Lba))
    return L_schur

def TB_to_harmonic_extension(T, B, leaves, internal):
    nleaves = len(leaves)
    Lbb = TB_to_L_block(T, B, internal, internal)
    Lba = TB_to_L_block(T, B, internal, leaves)
    L_schur = TB_to_L_schur(T, B, leaves)
    w, v1 = scipy.linalg.eigh(L_schur, eigvals=(1, nleaves-1))
    v2 = -np.dot(np.dot(np.linalg.pinv(Lbb), Lba), v1)
    v = np.vstack([v1, v2])
    return w, v

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

def T_to_outside_in(T):
    """
    Get an ordered sequence of vertices.
    This is according to a certain notion of centrality.
    The leaves are guaranteed to be first.
    @return: a list of ordered vertices
    """
    v_to_cent = T_to_v_to_centrality(T)
    pairs = sorted((c, v) for v, c in v_to_cent.items())
    return [v for c, v in pairs]

def T_to_inside_out(T):
    """
    Get an ordered sequence of vertices.
    This is according to a certain notion of centrality.
    The initial subsequences are guaranteed to induce connected graphs.
    The leaves are guaranteed to be last.
    @return: a list of ordered vertices
    """
    return list(reversed(T_to_outside_in(T)))

def T_to_order(T):
    """
    Get an ordered sequence of vertices.
    The order should have the property that every graph induced by
    an initial subsequence should be connected.
    @return: a list of ordered vertices
    """
    return R_to_preorder(T_to_R_canonical(T))

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


# Testing

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

g_example_b_T = set([
    mkedge(0,1),
    mkedge(0,2),
    mkedge(0,6),
    mkedge(2,3),
    mkedge(2,4),
    mkedge(4,5)])

g_example_b_B = {
        mkedge(0,1) : 1,
        mkedge(0,2) : 2,
        mkedge(0,6) : 2,
        mkedge(2,3) : 3,
        mkedge(2,4) : 3,
        mkedge(4,5) : 3}


class TestFtree(unittest.TestCase):

    def test_leaves_a(self):
        observed = T_to_leaves(g_example_T)
        expected = [1, 4, 5, 7]
        self.assertEqual(observed, expected)

    def test_leaves_b(self):
        observed = T_to_leaves(g_example_b_T)
        expected = [1, 3, 5, 6]
        self.assertEqual(observed, expected)

    def test_internal_vertices_a(self):
        observed = T_to_internal_vertices(g_example_T)
        expected = [2, 3, 6]
        self.assertEqual(observed, expected)

    def test_internal_vertices_b(self):
        observed = T_to_internal_vertices(g_example_b_T)
        expected = [0, 2, 4]
        self.assertEqual(observed, expected)

    def test_distance_matrix(self):
        vertices = (3,4,5)
        observed = TB_to_D(g_example_T, g_example_B, vertices)
        expected = np.array([
            [0, 4, 3],
            [4, 0, 7],
            [3, 7, 0]], dtype=float)
        self.assertTrue(np.allclose(observed, expected))

    def test_tree_centrality(self):
        observed = T_to_outside_in(g_example_T)
        expected = [1, 4, 5, 7, 2, 6, 3]
        self.assertEqual(observed, expected)

    def test_incidence_decomposition(self):
        vertices = T_to_outside_in(g_example_T)
        B = T_to_incidence_matrix(g_example_T, vertices)
        W = TB_to_weight_matrix(g_example_T, g_example_B)
        observed = np.dot(B, np.dot(W, B.T))
        expected = TB_to_L_principal(g_example_T, g_example_B, vertices)
        self.assertTrue(np.allclose(observed, expected))

    def test_schur_to_distance(self):
        leaves = T_to_leaves(g_example_T)
        # Compute the Schur complement Laplacian and the leaf distance matrix.
        L_schur = TB_to_L_schur(g_example_T, g_example_B, leaves)
        Dpp_direct = TB_to_D(g_example_T, g_example_B, leaves)
        # Compute one from the other.
        HDppH_schur = -2*np.linalg.pinv(L_schur)
        d = np.diag(HDppH_schur)
        e = np.ones_like(d)
        Dpp_schur = HDppH_schur - 0.5*(np.outer(d,e) + np.outer(e,d))
        self.assertTrue(np.allclose(Dpp_schur, Dpp_direct))

    def test_distance_to_schur(self):
        leaves = T_to_leaves(g_example_T)
        # Compute the Schur complement Laplacian and the leaf distance matrix.
        L_schur = TB_to_L_schur(g_example_T, g_example_B, leaves)
        Dpp_direct = TB_to_D(g_example_T, g_example_B, leaves)
        # Compute one from the other.
        HDppH = MatrixUtil.double_centered(Dpp_direct)
        L_schur_estimate = np.linalg.pinv(-0.5*HDppH)
        self.assertTrue(np.allclose(L_schur_estimate, L_schur))

    def test_L_schur_spectrum(self):
        leaves_a = T_to_leaves(g_example_T)
        leaves_b = T_to_leaves(g_example_b_T)
        L_schur_a = TB_to_L_schur(g_example_T, g_example_B, leaves_a)
        L_schur_b = TB_to_L_schur(g_example_b_T, g_example_b_B, leaves_b)
        w_a = scipy.linalg.eigh(L_schur_a, eigvals_only=True)
        w_b = scipy.linalg.eigh(L_schur_b, eigvals_only=True)
        self.assertTrue(np.allclose(w_a, w_b))

    def test_harmonic_extension_a(self):
        leaves = T_to_leaves(g_example_T)
        internal = T_to_internal_vertices(g_example_T)
        w_observed, v_observed = TB_to_harmonic_extension(
                g_example_T, g_example_B, leaves, internal)
        w_expected = [0.1707228213, 0.271592036629, 0.684669269055]
        self.assertTrue(np.allclose(w_observed, w_expected))

    def test_harmonic_extension_b(self):
        leaves = T_to_leaves(g_example_b_T)
        internal = T_to_internal_vertices(g_example_b_T)
        w_observed, v_observed = TB_to_harmonic_extension(
                g_example_b_T, g_example_b_B, leaves, internal)
        w_expected = [0.1707228213, 0.271592036629, 0.684669269055]
        self.assertTrue(np.allclose(w_observed, w_expected))

    def test_harmonic_extension_c(self):
        T = set([
            frozenset([2, 4]), frozenset([0, 6]), frozenset([2, 3]),
            frozenset([0, 2]), frozenset([0, 1]), frozenset([4, 5])])
        B = {
                frozenset([2, 4]): 3.0, frozenset([0, 6]): 2.0,
                frozenset([2, 3]): 3.0, frozenset([0, 2]): 2.0,
                frozenset([0, 1]): 1.0, frozenset([4, 5]): 3.0}
        leaves = T_to_leaves(T)
        internal = T_to_internal_vertices(T)
        w_observed, v_observed = TB_to_harmonic_extension(
                T, B, leaves, internal)
        w_expected = [0.1707228213, 0.271592036629, 0.684669269055]
        self.assertTrue(np.allclose(w_observed, w_expected))

    def test_R_to_vertex_partition(self):
        R = set([(8, 5), (8, 6), (8, 7), (6, 1), (6, 2), (7, 3), (7, 4)])
        observed = R_to_vertex_partition(R)
        expected = set([frozenset([5]), frozenset([1,2,6]), frozenset([3,4,7])])
        self.assertEqual(observed, expected)


if __name__ == '__main__':
    unittest.main()

