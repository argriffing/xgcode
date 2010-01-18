"""
Stuff related to spanning trees in various dimensions on various surfaces.
"""

import unittest
import heapq

def kruskal(V, E):
    """
    Pseudocode from wikipedia Kruskal's algorithm:
    1   function Kruskal(G)
    2     for each vertex v in G do
    3       Define an elementary cluster C(v) <-- {v}.
    4     Initialize a priority queue Q to contain all edges in G, using the weights as keys.
    5     Define a tree T <-- empty set       //T will ultimately contain the edges of the MST
    6     // n is total number of vertices
    7     while T has fewer than n-1 edges do
    8       // edge u,v is the minimum weighted route from/to v
    9       (u,v) <-- Q.removeMin()
    10      // prevent cycles in T. add u,v only if T does not already contain an edge consisting of u and v. 
    11      // Note that the cluster contains more than one vertex only if an edge containing a pair of
    12      // the vertices has been added to the tree.
    13      Let C(v) be the cluster containing v, and let C(u) be the cluster containing u.
    14      if C(v) != C(u) then
    15        Add edge (v,u) to T.
    16        Merge C(v) and C(u) into one cluster, that is, union C(v) and C(u).
    17    return tree T
    @param V: a list of hashable vertices
    @param E: a list of hashable (nonnegative weight, vertex u, vertex v) triples
    """
    # validate the input
    for weight, a, b in E:
        if (a not in V) or (b not in V):
            raise ValueError('expected both endpoints of an edge to be valid vertices')
    # define an elementary cluster for each vertex
    C = dict((v, set([v])) for v in V)
    # initialize a priority queue Q to contain all edges in G, using the weights as keys.
    Q = E[:]
    heapq.heapify(Q)
    # define a tree that will ultimately contain the edges of the minimum spanning tree
    T = set()
    # n is the total number of vertices
    n = len(V)
    while len(T) < n-1:
        # edge (u, v) is the minimum weighted route from/to v
        edge = heapq.heappop(Q)
        weight, u, v = edge
        if C[v] != C[u]:
            T.add(edge)
            C[u].update(C[v])
            for vertex in C[u]:
                C[vertex] = C[u]
    return T


class TestMST(unittest.TestCase):

    def test_kruskal(self):
        """
        An example of Kruskal's algorithm on Wikipedia.
        """
        V = set('abcdefg')
        E = [
            (7, 'a', 'b'),
            (5, 'a', 'd'),
            (9, 'd', 'b'),
            (8, 'b', 'c'),
            (7, 'b', 'e'),
            (15, 'd', 'e'),
            (6, 'd', 'f'),
            (5, 'c', 'e'),
            (9, 'e', 'g'),
            (11, 'f', 'g'),
            (8, 'f', 'e')
            ]
        T = kruskal(V, E)
        E_MST = set([
            (7, 'a', 'b'),
            (5, 'a', 'd'),
            (7, 'b', 'e'),
            (6, 'd', 'f'),
            (5, 'c', 'e'),
            (9, 'e', 'g')
            ])
        self.assertEqual(set(T), set(E_MST))


if __name__ == '__main__':
    unittest.main()
