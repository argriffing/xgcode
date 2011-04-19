"""
This is a function oriented rather than object oriented approach to trees.
"""

import unittest
from StringIO import StringIO
from collections import defaultdict

import NewickIO
import Ftree


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

def R_to_newick(R):
    """
    @param R: a directed topology
    @return: a newick string
    """
    r = Ftree.R_to_root(R)
    return _v_to_newick(Ftree.get_v_to_sinks(R), r) + ';'

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
        blen = B[Ftree.mkedge(v, v_to_source[v])]
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
    r = Ftree.R_to_root(R)
    v_to_source = Ftree.get_v_to_source(R)
    v_to_sinks = Ftree.get_v_to_sinks(R)
    return _Bv_to_newick(v_to_source, v_to_sinks, B, r) + ';'

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
    return Ftree.R_to_T(R)

def newick_to_TB(s, name_type=None):
    """
    @param s: newick string
    @return: undirected topology, branch lengths
    """
    R, B = newick_to_RB(s, name_type)
    return Ftree.R_to_T(R), B

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
    Ftree.mkedge(2,1),
    Ftree.mkedge(2,3),
    Ftree.mkedge(2,4),
    Ftree.mkedge(3,5),
    Ftree.mkedge(3,6),
    Ftree.mkedge(6,7)])

g_example_B = {
        Ftree.mkedge(2,1) : 1,
        Ftree.mkedge(2,3) : 2,
        Ftree.mkedge(2,4) : 2,
        Ftree.mkedge(3,5) : 3,
        Ftree.mkedge(3,6) : 3,
        Ftree.mkedge(6,7) : 3}

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

if __name__ == '__main__':
    unittest.main()

