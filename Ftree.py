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

import NewickIO
import MatrixUtil

def invseq(seq):
    """
    @param seq: an iterable
    @return: a map from the value to the sequence index
    """
    return dict((v, i) for i, v in enumerate(seq))

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

def R_to_T(R):
    """
    @param R: a directed topology
    @return: an undirected topology
    """
    return set(frozenset(d_edge) for d_edge in R)

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


# Newick IO stuff

class _IO_Node:
    """
    Force this to work with the NewickIO parser.
    """
    def __init__(self, name_type):
        self.name_type = name_type
        self.name = None
        self.blen = None
        self.children = []
    def __call__(self):
        return _IO_Node(self.name_type)
    def set_parent(self, node):
        pass
    def set_branch_length(self, blen):
        self.blen = blen
    def add_child(self, child):
        self.children.append(child)
    def add_name(self, name):
        if self.name_type:
            self.name = self.name_type(name)
        else:
            self.name = name

class _IO_Tree:
    """
    Force this to work with the NewickIO parser.
    The functions prefixed with x_ are not part
    of the NewickIO parser interface.
    """
    def __init__(self, name_type):
        self.root = None
        self.NodeFactory = _IO_Node(name_type)
        self.name_type = name_type
    def __call__(self):
        return _IO_Tree(self.name_type)
    def set_root(self, node):
        self.root = node
    def x_all_nodes(self):
        nodes = []
        shell = [self.root]
        while shell:
            nodes.extend(shell)
            next_shell = []
            for node in shell:
                next_shell.extend(node.children)
            shell = next_shell
        return nodes
    def x_get_RB(self):
        R = set()
        B = {}
        shell = [self.root]
        while shell:
            next_shell = []
            for node in shell:
                for child in node.children:
                    d_edge = (node.name, child.name)
                    R.add(d_edge)
                    B[frozenset(d_edge)] = child.blen
                    next_shell.append(child)
            shell = next_shell
        return R, B
    def x_assert_root(self):
        if self.root is None:
            raise ValueError('no root was specified')
    def x_assert_names(self):
        nodes = self.x_all_nodes()
        if any(n.name is None for n in nodes):
            msg = 'all nodes should be named, including internal nodes'
            raise ValueError(msg)
        if len(set(n.name for n in nodes)) < len(nodes):
            msg = 'nodes should have unique names'
            raise ValueError(msg)
    def x_assert_branch_lengths(self):
        root = self.root
        if root.blen is not None:
            raise ValueError('the root should not have a branch length')
        nodes = self.x_all_nodes()
        non_root_nodes = [n for n in nodes if n is not root]
        if any(n.blen is None for n in non_root_nodes):
            msg = 'every node except the root should have a branch length'
            raise ValueError(msg)

def newick_to_T(s, name_type=None):
    """
    Everything to do with branch lengths is ignored.
    @param s: newick string
    @return: undirected topology
    """
    tree = NewickIO.parse(s, _IO_Tree(name_type))
    tree.x_assert_root()
    tree.x_assert_names()
    R, B = tree.x_get_RB()
    return R_to_T(R)

def newick_to_TB(s, name_type=None):
    """
    @param s: newick string
    @return: undirected topology, branch lengths
    """
    R, B = newick_to_RB(s, name_type)
    return R_to_T(R), B

def newick_to_RB(s, name_type=None):
    """
    @param s: newick string
    @return: directed topology, branch lengths
    """
    tree = NewickIO.parse(s, _IO_Tree(name_type))
    tree.x_assert_root()
    tree.x_assert_names()
    tree.x_assert_branch_lengths()
    return tree.x_get_RB()


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

    def test_unweighted_from_newick(self):
        s = '(1, (5, (7)6)3, 4)2;'
        observed = newick_to_T(s, int)
        expected = g_example_T
        self.assertEqual(observed, expected)

    def test_weighted_from_newick(self):
        s = '(1:1, (5:3, (7:3)6:3)3:2, 4:2)2;'
        observed_T, observed_B = newick_to_TB(s, int)
        expected_T = g_example_T
        expected_B = g_example_B
        self.assertEqual(observed_T, expected_T)
        self.assertEqual(observed_B, expected_B)

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
        internal = T_to_internal_vertices(g_example_T)
        Lpp = TB_to_L_block(g_example_T, g_example_B, leaves, leaves)
        Lpr = TB_to_L_block(g_example_T, g_example_B, leaves, internal)
        Lrp = TB_to_L_block(g_example_T, g_example_B, internal, leaves)
        Lrr = TB_to_L_block(g_example_T, g_example_B, internal, internal)
        # Compute the Schur complement Laplacian and the leaf distance matrix.
        L_schur = Lpp - np.dot(Lpr, np.dot(np.linalg.pinv(Lrr), Lrp))
        Dpp_direct = TB_to_D(g_example_T, g_example_B, leaves)
        # Compute one from the other.
        HDppH_schur = -2*np.linalg.pinv(L_schur)
        d = np.diag(HDppH_schur)
        e = np.ones_like(d)
        Dpp_schur = HDppH_schur - 0.5*(np.outer(d,e) + np.outer(e,d))
        self.assertTrue(np.allclose(Dpp_schur, Dpp_direct))

    def test_distance_to_schur(self):
        leaves = T_to_leaves(g_example_T)
        internal = T_to_internal_vertices(g_example_T)
        Lpp = TB_to_L_block(g_example_T, g_example_B, leaves, leaves)
        Lpr = TB_to_L_block(g_example_T, g_example_B, leaves, internal)
        Lrp = TB_to_L_block(g_example_T, g_example_B, internal, leaves)
        Lrr = TB_to_L_block(g_example_T, g_example_B, internal, internal)
        # Compute the Schur complement Laplacian and the leaf distance matrix.
        L_schur = Lpp - np.dot(Lpr, np.dot(np.linalg.pinv(Lrr), Lrp))
        Dpp_direct = TB_to_D(g_example_T, g_example_B, leaves)
        # Compute one from the other.
        HDppH = MatrixUtil.double_centered(Dpp_direct)
        L_schur_estimate = np.linalg.pinv(-0.5*HDppH)
        self.assertTrue(np.allclose(L_schur_estimate, L_schur))


if __name__ == '__main__':
    unittest.main()

