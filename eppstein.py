"""
Implement parts of a tree layout algorithm described by Carlson and Eppstein.
"""

from collections import defaultdict
import unittest
import math

import graph
import petree

SUBTREE_PATH = 'path-subtree'
SUBTREE_RAKE = 'rake-subtree'
SUBTREE_OTHER = 'complicated-subtree'

TREE_TYPE_PATH = 'path-tree'
TREE_TYPE_RAKE = 'rake-tree'
TREE_TYPE_TRIPLE_RAKE = 'triple-rake-tree'
TREE_TYPE_OTHER = 'complicated-tree'

def get_plane_embedded_triple_rake_ar(s, d):
    """
    Get the angular resolution.
    @param s: the number of short paths
    @param d: the number of double turns
    """
    return math.pi * (0.5 + 1.0 / (2.0*(9 - 2*s + 2*d)))

def get_dp_edges(V, U):
    """
    Return a sequence of directed edges of an unrooted tree.
    The returned sequence of directed edges is sorted for dynamic programming.
    When a directed edge appears in the list,
    all of the directed edges it points towards will have preceded it.
    @param V: vertices
    @param U: undirected edges
    @return: a sequence of directed edges
    """
    V_meta = set()
    D_meta = set()
    nd = graph.g_to_nd(V, U)
    for v, neighbors in nd.items():
        for a in neighbors:
            V_meta.add((a, v))
            V_meta.add((v, a))
            for b in neighbors:
                if a != b:
                    D_meta.add(((a, v), (v, b)))
    return graph.topo_sort(V_meta, D_meta)

def get_tree_type(V, U):
    max_degree = graph.g_to_max_degree(V, U)
    if max_degree < 3:
        return TREE_TYPE_PATH
    elif max_degree == 3:
        nd = graph.g_to_nd(V, U)
        edge_to_max_degree = get_subtree_max_degree_map(V, U)
        # Count the number of vertices
        # that see a degree three vertex in three directions.
        ncenters = 0
        for v, neighbors in nd.items():
            nvisible = 0
            for w in neighbors:
                if edge_to_max_degree[(v, w)] == 3:
                    nvisible += 1
            if nvisible == 3:
                ncenters += 1
        if ncenters == 0:
            return TREE_TYPE_RAKE
        elif ncenters == 1:
            return TREE_TYPE_TRIPLE_RAKE
        else:
            return TREE_TYPE_OTHER
    else:
        return TREE_TYPE_OTHER

def get_high_degree_mst(V, U):
    """
    Get the high degree vertices and the paths connecting them.
    In this context high degree means at least three.
    Every high degree vertex is in this minimum spanning tree.
    Every vertex that sees at least two high degree vertices
    is in this minimum spanning tree.
    @return: V_mst, U_mst
    """
    V_mst = set()
    U_mst = set()
    nd = graph.g_to_nd(V, U)
    edge_to_max_degree = get_subtree_max_degree_map(V, U)
    for v, neighbors in nd.items():
        if len(nd[v]) > 2:
            V_mst.add(v)
        nvisible = 0
        for w in neighbors:
            if edge_to_max_degree[(v, w)] > 2:
                nvisible += 1
        if nvisible > 1:
            V_mst.add(v)
    U_mst = set(u for u in U if u <= V_mst)
    return V_mst, U_mst



def get_subtree_max_degree_map(V, U):
    """
    @param V: vertices
    @param U: undirected edges
    @return: a map from a directed edge to a max vertex degree
    """
    nd = graph.g_to_nd(V, U)
    dp_edges = get_dp_edges(V, U)
    edge_to_max_degree = dict()
    for a, b in reversed(dp_edges):
        neighbors = nd[b]
        degrees = set([len(neighbors)])
        for c in neighbors:
            if c != a:
                degrees.add(edge_to_max_degree[(b, c)])
        edge_to_max_degree[(a, b)] = max(degrees)
    return edge_to_max_degree

def get_subtree_type_map(V, U):
    """
    @param V: vertices
    @param U: undirected edges
    @return: a map from a directed edge to a subtree type
    """
    nd = graph.g_to_nd(V, U)
    dp_edges = get_dp_edges(V, U)
    edge_to_type = dict()
    for a, b in reversed(dp_edges):
        neighbors = nd[b]
        if len(neighbors) == 1:
            # an edge leading to a leaf is a path
            edge_to_type[(a, b)] = SUBTREE_PATH
        elif len(neighbors) == 2:
            # an extra edge at the root of a subtree does not change the type
            for c in neighbors:
                if c != a:
                    edge_to_type[(a, b)] = edge_to_type[(b, c)]
        elif len(neighbors) == 3:
            # trifurcations are the interesting case
            typecount = defaultdict(int)
            for c in neighbors:
                if c != a:
                    typecount[edge_to_type[(b, c)]] += 1
            if SUBTREE_OTHER in typecount:
                t = SUBTREE_OTHER
            elif SUBTREE_RAKE in typecount:
                if typecount[SUBTREE_RAKE] > 1:
                    t = SUBTREE_OTHER
                else:
                    t = SUBTREE_RAKE
            elif typecount[SUBTREE_PATH] > 1:
                t = SUBTREE_RAKE
            else:
                t = SUBTREE_PATH
            edge_to_type[(a, b)] = t
        else:
            # multifurcations of high degree imply neither path nor rake
            edge_to_type[(a, b)] = SUBTREE_OTHER
    return edge_to_type

def get_triple_rake_short_path_count(V, U):
    """
    @return: the number of short paths either 0, 1, 2, or 3
    """
    nd = graph.g_to_nd(V, U)
    edge_to_max_degree = get_subtree_max_degree_map(V, U)
    # Look for the vertex that sees a degree three vertex in three directions.
    vertex_to_nvisible = dict((v, 0) for v in V)
    for v, neighbors in nd.items():
        for w in neighbors:
            if edge_to_max_degree[(v, w)] == 3:
                vertex_to_nvisible[v] += 1
        if vertex_to_nvisible[v] == 3:
            center = v
    # Check the three directions for short paths.
        
        


class TestEppstein(unittest.TestCase):

    def test_subtree_type(self):
        V = set(range(7))
        U = set([
            graph.fs(5, 0),
            graph.fs(5, 1),
            graph.fs(5, 4),
            graph.fs(4, 2),
            graph.fs(4, 3),
            graph.fs(4, 6)])
        observed = get_subtree_type_map(V, U)
        expected = {
                (6, 4): SUBTREE_OTHER,
                (5, 4): SUBTREE_OTHER,
                (4, 3): SUBTREE_PATH,
                (4, 5): SUBTREE_RAKE,
                (1, 5): SUBTREE_OTHER,
                (0, 5): SUBTREE_OTHER,
                (5, 0): SUBTREE_PATH,
                (5, 1): SUBTREE_PATH,
                (4, 2): SUBTREE_PATH,
                (3, 4): SUBTREE_OTHER,
                (2, 4): SUBTREE_OTHER,
                (4, 6): SUBTREE_PATH}
        self.assertEqual(observed, expected)

    def test_subtree_max_degree(self):
        V = set(range(7))
        U = set([
            graph.fs(5, 0),
            graph.fs(5, 1),
            graph.fs(5, 4),
            graph.fs(4, 2),
            graph.fs(4, 3),
            graph.fs(4, 6)])
        observed = get_subtree_max_degree_map(V, U)
        expected = {
                (6, 4): 4,
                (5, 4): 4,
                (4, 3): 1,
                (4, 5): 3,
                (1, 5): 4,
                (0, 5): 4,
                (5, 0): 1,
                (5, 1): 1,
                (4, 2): 1,
                (3, 4): 4,
                (2, 4): 4,
                (4, 6): 1}
        self.assertEqual(observed, expected)

    def test_tree_type_high_degree(self):
        V = set(range(7))
        U = set([
            graph.fs(5, 0),
            graph.fs(5, 1),
            graph.fs(5, 4),
            graph.fs(4, 2),
            graph.fs(4, 3),
            graph.fs(4, 6)])
        observed = get_tree_type(V, U)
        self.assertEqual(observed, TREE_TYPE_OTHER)

    def test_tree_type_rake(self):
        V = set(range(6))
        U = set([
            graph.fs(5, 0),
            graph.fs(5, 1),
            graph.fs(5, 4),
            graph.fs(4, 2),
            graph.fs(4, 3)])
        observed = get_tree_type(V, U)
        self.assertEqual(observed, TREE_TYPE_RAKE)

    def test_tree_type_path(self):
        V = set(range(6))
        U = set([
            graph.fs(5, 0),
            graph.fs(0, 3),
            graph.fs(3, 2),
            graph.fs(1, 4),
            graph.fs(2, 1)])
        observed = get_tree_type(V, U)
        self.assertEqual(observed, TREE_TYPE_PATH)

    def test_tree_type_triple_rake(self):
        V = set(range(18))
        U = set([
            graph.fs(0, 1),
            graph.fs(1, 2),
            graph.fs(2, 3),
            graph.fs(3, 4),
            graph.fs(4, 5),
            graph.fs(5, 6),
            graph.fs(6, 7),
            graph.fs(7, 8),
            #
            graph.fs(7, 9),
            graph.fs(6, 10),
            graph.fs(5, 11),
            graph.fs(3, 12),
            graph.fs(2, 13),
            graph.fs(1, 14),
            #
            graph.fs(4, 15),
            graph.fs(15, 16),
            graph.fs(15, 17)])
        observed = get_tree_type(V, U)
        self.assertEqual(observed, TREE_TYPE_TRIPLE_RAKE)

    def test_tree_type_other(self):
        V = set(range(18))
        U = set([
            graph.fs(0, 1),
            #
            #graph.fs(1, 2),
            graph.fs(1, 10),
            #
            graph.fs(2, 3),
            graph.fs(3, 4),
            graph.fs(4, 5),
            graph.fs(5, 6),
            graph.fs(6, 7),
            graph.fs(7, 8),
            #
            graph.fs(7, 9),
            graph.fs(6, 10),
            graph.fs(5, 11),
            graph.fs(3, 12),
            graph.fs(2, 13),
            graph.fs(1, 14),
            #
            graph.fs(4, 15),
            graph.fs(15, 16),
            graph.fs(15, 17)])
        observed = get_tree_type(V, U)
        self.assertEqual(observed, TREE_TYPE_OTHER)

    def test_high_degree_mst(self):
        V = set(range(18)) - set([12])
        U = set([
            graph.fs(0, 1),
            graph.fs(1, 2),
            graph.fs(2, 3),
            graph.fs(3, 4),
            graph.fs(4, 5),
            graph.fs(5, 6),
            graph.fs(6, 7),
            graph.fs(7, 8),
            #
            graph.fs(7, 9),
            graph.fs(6, 10),
            graph.fs(5, 11),
            #
            #graph.fs(3, 12),
            #
            graph.fs(2, 13),
            graph.fs(1, 14),
            #
            graph.fs(4, 15),
            graph.fs(15, 16),
            graph.fs(15, 17)])
        V_mst, U_mst = get_high_degree_mst(V, U)
        V_mst_expected = set([1, 2, 3, 4, 5, 6, 7, 15])
        U_mst_expected = set([
            graph.fs(1, 2),
            graph.fs(2, 3),
            graph.fs(3, 4),
            graph.fs(4, 5),
            graph.fs(5, 6),
            graph.fs(6, 7),
            graph.fs(4, 15)])
        self.assertEqual(V_mst, V_mst_expected)
        self.assertEqual(U_mst, U_mst_expected)


if __name__ == '__main__':
    unittest.main()

