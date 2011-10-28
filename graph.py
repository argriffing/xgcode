"""
This is a generic graph module.

I am audacious enough to think that this will be better
than NetworkX for my purposes.
Like NetworkX, the core objects are distinct hashable vertices.
Unlike NetworkX, the graph is not implemented through dictionaries.
Some notational convention follows.
V is the set of (distinct hashable) vertices in the graph.
E is the set of undirected edges,
where an unordered edge is a frozenset of two distinct vertices.
D is the set of directed edges,
where a directed edge is an (ordered) tuple of two distinct vertices.
nd is the neighbor dictionary that maps a vertex its set of neighbors.
pd is the parent dictionary that maps a vertex to its set of parents.
cd is the child dictionary that maps a vertex to its set of children.
The following notation is used in function names.
g is an undirected graph.
dg is a directed graph.
dag is a directed acyclic graph.
fs is frozenset.
"""

import unittest

class TopoError(Exception): pass

def fs(*args):
    return frozenset(args)

def topo_sort(V, D):
    """
    Topologically sort a directed acyclic graph.
    Use an algorithm from http://en.wikipedia.org/wiki/Topological_sorting
    @param V: set of vertices
    @param D_in: set of directed edges
    @return: topologically sorted sequence of vertices
    """
    cd = dag_to_cd(V, D)
    pd = dag_to_pd(V, D)
    L = []
    S = set(v for v in V if not pd[v])
    while S:
        n = S.pop()
        L.append(n)
        children = set(cd[n])
        for m in children:
            # remove the edge from the graph
            cd[n].remove(m)
            pd[m].remove(n)
            # check for other incoming edges to m
            if not pd[m]:
                S.add(m)
    if any(cd.values()):
        raise TopoError('the graph has a directed cycle')
    return L

def g_to_max_degree(V, E):
    """
    Given an undirected graph compute its max degree.
    """
    nd = g_to_nd(V, E)
    return max(len(neighbors) for neighbors in nd.values())

def g_to_nd(V, E):
    """
    Given an undirected graph compute the neighbor dictionary.
    A vertex will exist as a key in the dictionary even if it has no neighbors.
    """
    nd = dict((v, set()) for v in V)
    for a, b in E:
        nd[a].add(b)
        nd[b].add(a)
    return nd

def nd_to_dag_component(nd, source):
    """
    Pick a single connected component out of an undirected graph.
    @param nd: neighbor dictionary for the whole graph
    @param source: a single vertex in the connected component of interest
    @return: V_component, D_component
    """
    D = set()
    shell = set([source])
    visited = set(shell)
    while shell:
        for a in shell:
            nshell = set()
            for b in nd[a]:
                if b not in visited:
                    D.add((a, b))
                    nshell.add(b)
                    visited.add(b)
            shell = nshell
    return visited, D

def g_to_dag_component(V, E, source):
    """
    Pick a single connected component out of an undirected graph.
    @param V: vertex set
    @param E: undirected edges
    @return: V_component, D_component
    """
    return nd_to_dag_component(g_to_nd(V, E), source)

def dag_to_cd(V, D):
    """
    @param V: vertex set
    @param D: directed edges
    @return: dictionary mapping a vertex to the set of child vertices
    """
    cd = dict((v, set()) for v in V)
    for a, b in D:
        cd[a].add(b)
    return cd

def dag_to_pd(V, D):
    """
    @param V: vertex set
    @param D: directed edges
    @return: dictionary mapping a vertex to the set of parent vertices
    """
    pd = dict((v, set()) for v in V)
    for a, b in D:
        pd[b].add(a)
    return pd



class TestGraph(unittest.TestCase):

    def test_graph(self):
        D = set([
            (5, 0), (5, 1), (5, 4),
            (4, 2), (4, 3)])
        V = set(range(6))
        E = set(fs(*p) for p in D)
        V_observed, D_observed = nd_to_dag_component(g_to_nd(V, E), 5)
        self.assertEqual(D, D_observed)

    def test_topo_sort_cycle(self):
        V = set([1,2,3])
        D = set([(1,2), (2,3), (3,1)])
        self.assertRaises(TopoError, topo_sort, V, D)

    def test_topo_sort_path(self):
        V = set([1,2,3])
        D = set([(1,3), (3,2)])
        self.assertEqual(topo_sort(V, D), [1, 3, 2])

    def test_topo_sort_fork(self):
        V = set([1,2,3,4])
        D = set([(4,3), (3,2), (3,1)])
        valid_sorts = [
                [4, 3, 2, 1],
                [4, 3, 1, 2]]
        self.assertIn(topo_sort(V, D), valid_sorts)

if __name__ == '__main__':
    unittest.main()

