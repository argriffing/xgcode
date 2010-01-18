"""Convert a tree to a directed acyclic graph on the directed edges of the tree. [UNFINISHED]
"""

import unittest

import Newick


class Node:
    """
    A node in a directed acyclic graph.
    """

    def __init__(self):
        """
        Initialize the forward and backward links.
        """
        self.next = []
        self.prev = []

    def link_to(self, next):
        """
        Adds a forward link to the target and a backward link from the target to this node.
        """
        self.next.append(next)
        next.prev.append(self)


class DAG:
    """
    A directed acyclic graph.
    A source is a node with no backward links.
    A sink is a node with no forward links.
    """

    def __init__(self):
        """
        Initialize the list of sources.
        """
        self.sources = []

    def add_source(self, source):
        """
        Append a source node to the list of sources.
        """
        self.sources.append(source)

    def preorder(self):
        """
        Yield nodes according to the partial order induced by the directed edges.
        """
        # Track the number of yielded dependencies for each node.
        # A node can be yielded after all of its dependencies have been yielded.
        id_to_count = {}
        # initialize the list of nodes that have no dependencies.
        sources = self.sources
        while sources:
            next_sources = []
            for node in sources:
                for target in node.next:
                    tid = id(target)
                    # count the number of dependencies accounted for in the target node
                    accounted_deps = id_to_count.get(tid, 0) + 1
                    # add one to the number of accounted-for dependencies in the target
                    id_to_count[tid] = accounted_deps
                    # if the target has all dependencies accounted for then it can be yielded
                    if accounted_deps == len(target.prev):
                        next_sources.append(target)
                    elif accounted_deps > len(target.prev):
                        assert False, 'maybe there is a cycle in the dependency structure'
                yield node
            sources = next_sources

    def postorder(self):
        return reversed(list(self.preorder()))


class TreeDAG(DAG):
    """
    Represent a tree by a directed acyclic graph on the directed edges of the tree.
    """

    def __init__(self, tree):
        """
        @param tree: a tree
        """
        pass


class TestTreeDAG(unittest.TestCase):

    def test_placeholder(self):
        pass


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestTreeDAG)
    unittest.TextTestRunner(verbosity=2).run(suite)

