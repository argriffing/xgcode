"""
This is a function oriented rather than object oriented approach to trees.

This is an interface module for the Ftree module.
In addition to the standard Ftree
variable names T, R, and B, we add N as the map from vertex to name.
With this new variable, vertices may now have names separate from
the hashable vertices themselves.
This facilitates, for example, the interpretation of a Newick string
whose leaves are labeled but whose internal vertices are not labeled.
"""

import unittest

#import dendropy

import NewickIO
import Ftree
from Ftree import mkedge

class FtreeIOError(Exception): pass

class _IO_Tree:
    """
    This implements the simplified NewickIO API.
    """
    def __init__(self):
        """
        Initialize the variables.
        """
        # This is the state of interest to the caller.
        self.v_to_name = {}
        self.B = {}
        self.R = set()
        # This is internal state.
        self.next_v = 0
        self.v_to_source = {}
        self.v_to_hanging_length = {}
    def create_root(self):
        root = self.next_v
        self.next_v += 1
        return root
    def add_child(self, parent, v):
        self.R.add((parent, v))
        self.v_to_source[v] = parent
        # resolve a hanging branch
        if v in self.v_to_hanging_length:
            self.set_branch_length(v, self.v_to_hanging_length[v])
            del self.v_to_hanging_length[v]
    def set_branch_length(self, v, blen):
        """
        Note that a branch length can be set to a root.
        This happens during the construction of the tree
        when the subtree has not yet been connected to the rest of the tree.
        """
        if v in self.v_to_source:
            edge = mkedge(v, self.v_to_source[v])
            self.B[edge] = blen
        else:
            self.v_to_hanging_length[v] = blen
    def set_name(self, v, name):
        self.v_to_name[v] = name
    def set_root(self, v):
        """
        This is slow, probably as a result of the design.
        """
        self.R = Ftree.T_to_R_specific(Ftree.R_to_T(self.R), v)
        self.v_to_source = Ftree.R_to_v_to_source(self.R)
    def finish(self):
        r = Ftree.R_to_root(self.R)
        if r in self.v_to_hanging_length:
            raise FtreeIOError('the root should not have a hanging branch')

def dendropy_to_RBN(tree):
    """
    @param tree: a dendropy.Tree object
    @return: R, B, N
    """
    next_v = 0
    node_to_v = {}
    R = set()
    B = {}
    N = {}
    for edge in tree.get_edge_set():
        # get the source and sink nodes and the branch length
        blen = edge.length
        head = edge.head_node
        tail = edge.tail_node
        if (head is None) or (tail is None):
            continue
        # add unknown vertices
        for node in (edge.head_node, edge.tail_node):
            if node not in node_to_v:
                v = next_v
                next_v += 1
                node_to_v[node] = v
                if node.taxon is not None:
                    N[v] = node.taxon.label
        # add the directed edge
        d_edge = (node_to_v[head], node_to_v[tail])
        R.add(d_edge)
        # add the branch length if it exists
        if blen is not None:
            B[d_edge] = blen
    return R, B, N

def R_to_newick(R):
    """
    @param R: a directed topology
    @return: a newick string
    """
    r = Ftree.R_to_root(R)
    return _v_to_newick(Ftree.R_to_v_to_sinks(R), r) + ';'

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

def _BNv_to_newick(v_to_source, v_to_sinks, B, N, v):
    """
    This is part of writing a newick string with branch lengths.
    Note that this function does not add the semicolon termination.
    @param v_to_source: a map from a vertex to its source
    @param v_to_sinks: a map from a vertex to a set of sinks
    @param B: branch lengths
    @param N: map from vertices to names
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
        if v not in N:
            raise ValueError('expected leaf vertices to be named')
        return str(N[v]) + suffix
    sinks = sorted(v_to_sinks[v])
    arr = [_BNv_to_newick(v_to_source, v_to_sinks, B, N, x) for x in sinks]
    internal_vertex_name = str(N.get(v, ''))
    return '(' + ', '.join(arr) + ')' + internal_vertex_name + suffix

def RB_to_newick(R, B):
    """
    @param R: a directed topology
    @param B: branch lengths
    @return: a newick string
    """
    r = Ftree.R_to_root(R)
    v_to_source = Ftree.R_to_v_to_source(R)
    v_to_sinks = Ftree.R_to_v_to_sinks(R)
    return _Bv_to_newick(v_to_source, v_to_sinks, B, r) + ';'

def RBN_to_newick(R, B, N):
    """
    @param R: a directed topology
    @param B: branch lengths
    @param N: map from vertices to names
    @return: a newick string
    """
    r = Ftree.R_to_root(R)
    v_to_source = Ftree.R_to_v_to_source(R)
    v_to_sinks = Ftree.R_to_v_to_sinks(R)
    return _BNv_to_newick(v_to_source, v_to_sinks, B, N, r) + ';'

def T_to_newick(T):
    """
    Get a newick string from an unweighted topology.
    @param T: topology
    @return: newick string
    """
    return R_to_newick(Ftree.T_to_R_canonical(T))

def TB_to_newick(T, B):
    """
    Get a newick string from a weighted topology.
    @param T: topology
    @param B: branch lengths
    @return: newick string
    """
    return RB_to_newick(Ftree.T_to_R_canonical(T), B)

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

def newick_to_TN(s):
    """
    Everything to do with branch lengths is ignored.
    @param s: newick string
    @return: undirected topology, vertex name map
    """
    tree = NewickIO.parse_simple(s, _IO_Tree())
    return Ftree.R_to_T(tree.R), tree.v_to_name

def newick_to_TBN(s):
    """
    @param s: newick string
    @return: undirected topology, branch lengths, vertex name map
    """
    tree = NewickIO.parse_simple(s, _IO_Tree())
    T = Ftree.R_to_T(tree.R)
    Ftree.TB_assert_branch_lengths(T, tree.B)
    return T, tree.B, tree.v_to_name

def newick_to_RBN(s):
    """
    @param s: newick string
    @return: directed topology, branch lengths, vertex name map
    """
    tree = NewickIO.parse_simple(s, _IO_Tree())
    Ftree.RB_assert_branch_lengths(tree.R, tree.B)
    return tree.R, tree.B, tree.v_to_name

def newick_to_T(s, name_type=None):
    """
    Everything to do with branch lengths is ignored.
    Vertex names are used as vertices.
    This is mostly for testing.
    @param s: newick string
    @return: undirected topology
    """
    T, N = newick_to_TN(s)
    N = get_validated_name_map(N, name_type)
    T = set(mkedge(N[a], N[b]) for a, b in T)
    return T

def newick_to_TB(s, name_type=None):
    """
    Vertex names are used as vertices.
    This is mostly for testing.
    @param s: newick string
    @return: undirected topology, branch lengths
    """
    T, B, N = newick_to_TBN(s)
    N = get_validated_name_map(N, name_type)
    T = set(mkedge(N[a], N[b]) for a, b in T)
    B = dict((mkedge(N[a], N[b]), x) for (a, b), x in B.items())
    return T, B

def newick_to_RB(s, name_type=None):
    """
    Vertex names are used as vertices.
    This is mostly for testing.
    @param s: newick string
    @return: directed topology, branch lengths
    """
    R, B, N = newick_to_RBN(s)
    N = get_validated_name_map(N, name_type)
    R = set((N[a], N[b]) for a, b in R)
    B = dict((mkedge(N[a], N[b]), x) for (a, b), x in B.items())
    return R, B

def newick_to_R(s, name_type=None):
    tree = NewickIO.parse_simple(s, _IO_Tree())
    N = get_validated_name_map(tree.v_to_name, name_type)
    return set((N[a], N[b]) for a, b in tree.R)

def get_validated_name_map(N, name_type):
    """
    Convert the vertex name map N to the requested type and check conditions.
    The conditions are that every vertex should have a name
    and that the names should be unique.
    @param N: vertex name map
    @param name_type: a function that defines the name type
    @return: a validated name map with names of the requested type
    """
    if name_type:
        N = dict((v, name_type(n)) for v, n in N.items())
    nvertices = len(N)
    names = N.values()
    if any(n is None for n in names):
        raise ValueError(
                'expected a name for each vertex, including internal vertices')
    if len(set(names)) < nvertices:
        raise ValueError('expected unique vertex names')
    return N


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

class TestFtreeIO(unittest.TestCase):

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

    def test_topo_b_from_newick(self):
        s = '((1:1, 2:0.5)6:1, (3:0.33333333333, 4:0.5)7:1, 5:1)8;'
        observed_T, observed_B = newick_to_TB(s, int)
        expected_T = Ftree.R_to_T(set([
            (8,7), (8,6), (8,5), (7,4), (7,3), (6,2), (6,1)]))
        self.assertEqual(observed_T, expected_T)
        observed_leaves = Ftree.T_to_leaves(observed_T)
        expected_leaves = [1, 2, 3, 4, 5]
        self.assertEqual(observed_leaves, expected_leaves)

    def _help_hanging_root_branch(self):
        """
        This is supposed to raise an exception.
        """
        s = '(1:1, (5:3, (7:3)6:3)3:2, 4:2)2:1;'
        T, B = newick_to_TB(s, int)

    def test_hanging_root_branch(self):
        self.assertRaises(FtreeIOError, self._help_hanging_root_branch)

    def test_dendropy_io(self):
        s = '((a:1, b:1), c:3, d:4)e;'
        tree = dendropy.Tree.get_from_string(s, 'newick')
        R, B, N = dendropy_to_RBN(tree)
        print R
        print B
        print N

if __name__ == '__main__':
    unittest.main()

