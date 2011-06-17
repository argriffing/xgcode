"""
This is a function oriented approach to sets of trees.

A forest will have a vertex set V, a set T of unordered edges,
and a dictionary B that maps unordered edges to branch lengths.
If matrices are used, they will be numpy arrays.
The initial motivation for this module is to partition leaves
according to the nodal domains associated with spectral graph theory.
A forest is basically a graph which may be connected
or unconnected and does not have cycles.
"""

import unittest
from collections import defaultdict


def mkedge(a, b):
    return frozenset([a, b])

def T_to_v_to_neighbors(T):
    """
    Note that if a vertex has no neighbors then it will not be in the map.
    @param T: set of undirected edges
    @return: a map from a vertex to the set of its neighboring vertices
    """
    d = defaultdict(set)
    for a, b in T:
        d[a].add(b)
        d[b].add(a)
    return d

def VT_to_vertex_partition(V, T):
    """
    @param V: set of vertices in the forest
    @param T: set of undirected edges
    """
    vertex_partition = set()
    # construct the vertex partition by getting tree vertex sets
    v_to_neighbors = T_to_v_to_neighbors(T)
    remaining = set(V)
    while remaining:
        shell = set([remaining.pop()])
        tree_vertices = set(shell)
        while shell:
            next_shell = set()
            for v in shell:
                for n in v_to_neighbors.get(v, ()):
                    if n not in tree_vertices:
                        next_shell.add(n)
            tree_vertices.update(next_shell)
            shell = next_shell
        vertex_partition.add(frozenset(tree_vertices))
    return vertex_partition


class TestFtree(unittest.TestCase):

    def test_vertex_partition(self):
        V = set([1,2,3,4,5,6,7])
        T = set([
            mkedge(2,3),
            mkedge(4,5),
            mkedge(4,6),
            mkedge(4,7),
            ])
        expected = set([
            frozenset([1]),
            frozenset([2,3]),
            frozenset([4,5,6,7]),
            ])
        observed = VT_to_vertex_partition(V, T)
        self.assertEqual(expected, observed)


if __name__ == '__main__':
    unittest.main()

