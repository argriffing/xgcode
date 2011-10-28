"""
This is a planar embedded tree module.

The planar embedded tree is associated with the following objects
a set of vertices,
a set of undirected edges,
a set of directed edges, and
a set of ordered pairs of directed edges.
The purpose of the extra planarity information
is for meaningful comparison among tree layout algorithms.
"""

import unittest

def _make_planar_tree_helper(
        labeled_trees, parent, parent_edge,
        V, U, D, P):
    """
    Make a planar embedded tree from plain input.
    @param labeled_trees: each is like (vertex, labeled_tree sequence)
    @param parent: the parent vertex
    @param parent_edge: None or the directed edge to the parent
    @param V: the set of integer vertices
    @param U: the set of undirected edges
    @param D: the set of directed edges
    @param P: the set of pairs of directed edges
    """
    v_prev = None
    for v, children in labeled_trees:
        V.add(v)
        edge = (parent, v)
        D.add(edge)
        D.add(tuple(reversed(edge)))
        U.add(frozenset(edge))
        _make_planar_tree_helper(
                children, v, (parent, v),
                V, U, D, P)
        if v_prev is not None:
            P.add(((v_prev, parent), (parent, v)))
        v_prev = v
    if parent_edge is not None:
        forward_parent_edge = parent_edge
        backward_parent_edge = tuple(reversed(parent_edge))
    if labeled_trees:
        forward_edge = (parent, labeled_trees[0][0])
        backward_edge = (labeled_trees[-1][0], parent)
        if parent_edge is not None:
            P.add((forward_parent_edge, forward_edge))
            P.add((backward_edge, backward_parent_edge))
        else:
            P.add((backward_edge, forward_edge))
    else:
        P.add((forward_parent_edge, backward_parent_edge))

def make_planar_tree(labeled_tree):
    """
    Make a planar embedded tree from an nx style rooted tree input.
    V is the set of integer vertices.
    U is the set of undirected edges.
    D is the set of directed edges.
    P is the set of pairs of directed edges.
    @param labeled_tree: like (vertex, labeled_tree sequence)
    @return: V, U, D, P
    """
    V = set()
    U = set()
    D = set()
    P = set()
    root, children = labeled_tree
    V.add(root)
    _make_planar_tree_helper(children, root, None, V, U, D, P)
    return V, U, D, P


class TestPlanarEmbeddedTree(unittest.TestCase):

    def test_petree(self):
        labeled_tree = [5, [
            [4, [
                [3, []],
                [2, []]]],
            [1, []],
            [0, []]]]
        D = set([
            (5, 0), (5, 1), (5, 4),
            (4, 2), (4, 3)])
        D.update(set(tuple(reversed(p)) for p in D))
        U = set(frozenset(p) for p in D)
        V = set(range(6))
        P = set([
            ((0,5), (5,4)),
            ((5,4), (4,3)),
            ((4,3), (3,4)),
            ((3,4), (4,2)),
            ((4,2), (2,4)),
            ((2,4), (4,5)),
            ((4,5), (5,1)),
            ((5,1), (1,5)),
            ((1,5), (5,0)),
            ((5,0), (0,5))])
        Vo, Uo, Do, Po = make_planar_tree(labeled_tree)
        self.assertEqual(V, Vo)
        self.assertEqual(U, Uo)
        self.assertEqual(D, Do)
        self.assertEqual(P, Po)


if __name__ == '__main__':
    unittest.main()

