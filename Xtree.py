"""
This module is supposed to conform to more mathematical definitions of X-tree variants.
"""

import unittest
import StringIO

import MatrixUtil

def edge_id(edge):
    """
    @param edge: a sequence of objects
    @return: a tuple of object ids
    """
    return tuple(id(element) for element in edge)


class EdgeInfo:
    """
    This is a simple utility aggregation of variables per directed edge.
    """

    def __init__(self):
        self.upstream_nknown = 0
        self.upstream_edges = []
        self.downstream_edges = []


class PXVertex:
    """
    This is the abstract base class of vertices of phylogenetic xtrees.
    It is for functions with the same implementation for both weighted and unweighted trees.
    """

    def has_label(self):
        return self.label is not None

    def get_neighbor(self):
        neighbors = self.get_neighbors()
        if len(neighbors) != 1:
            raise ValueError('expected a single neighbor but found %d' % len(neighbors))
        return neighbors[0]

    def get_exits(self, entrance):
        return [v for v in self.get_neighbors() if v is not entrance]

    def get_vertices(self):
        """
        @return: a list of vertices
        """
        return self.get_preorder_vertices()

    def get_preorder_vertices(self):
        """
        @return: a list of vertices in a partial order
        """
        vertices = [self]
        i = 0
        while i < len(vertices):
            vertices.extend(vertices[i].children)
            i += 1
        return vertices

    def get_labeled_vertices(self):
        """
        @return: a list of labeled vertices
        """
        return [vertex for vertex in self.get_vertices() if vertex.has_label()]

    def get_unlabeled_vertices(self):
        """
        @return: a list of unlabeled vertices
        """
        return [vertex for vertex in self.get_vertices() if not vertex.has_label()]

    def get_vertex_id_triples(self):
        """
        @return: a list of vertex id triples
        """
        triples = []
        for center in self.get_unlabeled_vertices():
            center_id = id(center)
            neighbor_ids = [id(v) for v in center.get_neighbors()]
            for a in neighbor_ids:
                for b in neighbor_ids:
                    if a != b:
                        triples.append((a, center_id, b))
        return triples

    def get_labels(self):
        """
        @return: the sorted list of labels
        """
        return list(sorted(v.label for v in self.get_labeled_vertices()))

    def get_nontrivial_splits(self):
        """
        Get label splits defined by the internal edges of the tree if self is the root.
        A split is a frozenset of two frozensets of labels.
        @return: a set of splits
        """
        return set(s for s in self.get_splits() if min(len(x) for x in s) > 1)

    def get_preorder_edges(self):
        """
        @return: an ordered list of ordered vertex pairs
        """
        preorder_edges = []
        # initialize the info per edge
        edge_to_info = dict((edge_id(edge), EdgeInfo()) for edge in self.get_unordered_edges())
        # define the dependence structure among the directed edges
        for a, b, c in self.get_vertex_id_triples():
            edge_to_info[(a,b)].downstream_edges.append((b,c))
            edge_to_info[(b,c)].upstream_edges.append((a,b))
        # initialize the label set for each directed edge with no upstream
        next_edges = []
        for v in self.get_labeled_vertices():
            edge = (v, v.get_neighbor())
            next_edges.append(edge_id(edge))
        # get the topological ordering of the directed acyclic graph
        while next_edges:
            preorder_edges.extend(next_edges)
            edges = next_edges
            next_edges = []
            for edge in edges:
                for downstream in edge_to_info[edge].downstream_edges:
                    ds_info = edge_to_info[downstream]
                    ds_info.upstream_nknown += 1
                    if ds_info.upstream_nknown == len(ds_info.upstream_edges):
                        next_edges.append(downstream)
        # return the vertex edges instead of the vertex id edges
        id_to_vertex = dict((id(v), v) for v in self.get_vertices())
        return [(id_to_vertex[a], id_to_vertex[b]) for a, b in preorder_edges]

    def get_postorder_edges(self):
        """
        @return: an ordered list of ordered vertex pairs
        """
        return list(reversed(self.get_preorder_edges()))

    def get_splits(self):
        """
        Get label splits defined by the edges of the tree if self is the root.
        A split is a frozenset of two frozensets of labels.
        @return: a set of splits
        """
        all_labels = set(self.get_labels())
        label_splits = set()
        edge_to_labels = {}
        for vertex_pair in self.get_postorder_edges():
            edge = edge_id(vertex_pair)
            source, sink = vertex_pair
            if sink.has_label():
                labels = set([sink.label])
            else:
                labels = set()
                for target in sink.get_exits(source):
                    # The labels of the next edge are known.
                    # This is because of the special ordering of the vertex pairs.
                    next_edge = edge_id((sink, target))
                    labels.update(edge_to_labels[next_edge])
            edge_to_labels[edge] = labels
            split = frozenset((frozenset(labels), frozenset(all_labels - labels)))
            label_splits.add(split)
        return label_splits

    def get_quartets(self):
        """
        Get the set of quartets defined by the tree.
        A quartet is a minimal nontrivial split with two labels on each side.
        The paper "An algorithm for reporting all shared quartets between two binary trees"
        by Thomas Mailund was helpful in implementing this algorithm.
        @return: a set of quartets
        """
        all_labels = set(self.get_labels())
        quartets = set()
        edge_to_labels = {}
        for vertex_pair in self.get_postorder_edges():
            edge = edge_id(vertex_pair)
            source, sink = vertex_pair
            if sink.has_label():
                labels = set([sink.label])
            else:
                labels = set()
                targets = sink.get_exits(source)
                # update the labels for dynamic programming
                for target in targets:
                    next_edge = edge_id((sink, target))
                    labels.update(edge_to_labels[next_edge])
                # get the quartets
                ordered_complement = list(all_labels - labels)
                for i, first_target in enumerate(targets):
                    for j, second_target in enumerate(targets):
                        if i < j:
                            first_edge = edge_id((sink, first_target))
                            second_edge = edge_id((sink, second_target))
                            for first_label in edge_to_labels[first_edge]:
                                for second_label in edge_to_labels[second_edge]:
                                    left_side = frozenset((first_label, second_label))
                                    for a, third_label in enumerate(ordered_complement):
                                        for b, fourth_label in enumerate(ordered_complement):
                                            if a < b:
                                                right_side = frozenset((third_label, fourth_label))
                                                quartet = frozenset((left_side, right_side))
                                                quartets.add(quartet)
            edge_to_labels[edge] = labels
        return quartets

    def get_newick_string(self):
        """
        @return: a newick string
        """
        return self.get_newick_substring() + ';'


class UPXVertex(PXVertex):
    """
    This is the vertex of an unweighted phylogenetic xtree.
    A phylogenetic xtree is defined mathematically in Semple and Steele.
    Only vertices of degree one are labeled.
    All vertices of degree one are labeled.
    Labels are unique.
    To keep the architecture light, this class should not know about branch lengths.
    """

    def __init__(self):
        self.parent = None
        self.children = []
        self.label = None

    def get_neighbors(self):
        parent_extension = [self.parent] if self.parent else []
        return self.children + parent_extension

    def add_child(self, child):
        self.children.append(child)
        child.parent = self

    def remove_child(self, child):
        child.parent = None
        self.children.remove(child)

    def reroot(self):
        """
        Reroot the tree at the current vertex.
        """
        # get the path to the root
        path = []
        v = self
        while v.parent:
            path.append((v.parent, v))
            v = v.parent
        # reverse parent child relationships on the path
        for parent, child in path:
            parent.remove_child(child)
            child.add_child(parent)

    def get_unordered_edges(self):
        """
        @return: a list of ordered vertex pairs
        """
        edges = []
        for v in self.get_vertices():
            if v.parent:
                edges.append((v, v.parent))
                edges.append((v.parent, v))
        return edges

    def get_distance_matrix(self):
        """
        This assumes that the labels are integers in [0, N).
        @return: a row major matrix
        """
        labels = self.get_labels()
        n = len(labels)
        if labels != range(n):
            raise ValueError()
        ordered_tips = [None]*n
        for vertex in self.get_labeled_vertices():
            ordered_tips[vertex.label] = vertex
        D = [[0]*n for i in range(n)]
        for i, vertex in enumerate(ordered_tips):
            distance = 0
            stack = [(vertex, vertex.get_neighbor(), 1)]
            while stack:
                source, sink, distance = stack.pop()
                if sink.has_label():
                    D[i][sink.label] = distance
                for next_sink in sink.get_neighbors():
                    if next_sink is not source:
                        stack.append((sink, next_sink, distance+1))
        return D

    def get_newick_substring(self):
        """
        @return: part of a newick string corresponding to the subtree
        """
        if self.has_label():
            return 'n%d' % self.label
        else:
            return '(' + ', '.join(child.get_newick_substring() for child in self.children) + ')'


class Branch:
    """
    This concept applies to weighted phylogenetic X-trees.
    To keep the architecture light, only the parent is referenced by the branch.
    This aggregation enforces the invariant that a vertex has a parent if and only if it has a branch length.
    """

    def __init__(self, parent, length):
        """
        @param parent: the parent vertex of the branch
        @param length: the length of the branch
        """
        self.parent = parent
        self.length = length


class WPXVertex(PXVertex):
    """
    A vertex of a weighted phylogenetic X-tree.
    """

    def __init__(self):
        self.branch = None
        self.children = []
        self.label = None

    def get_neighbors(self):
        parent_extension = [self.branch.parent] if self.branch else []
        return self.children + parent_extension

    def get_neighbors_and_branch_lengths(self):
        """
        @return: a list of (neighbor, branch_length) pairs
        """
        parent_extension = [(self.branch.parent, self.branch.length)] if self.branch else []
        return [(child, child.branch.length) for child in self.children] + parent_extension

    def get_neighbor_and_branch_length(self):
        """
        @return: a (neighbor, branch_length) pair
        """
        pairs = self.get_neighbors_and_branch_lengths()
        if len(pairs) != 1:
            raise ValueError('expected a single neighbor but found %d' % len(pairs))
        return pairs[0]

    def has_label(self):
        return self.label is not None

    def add_child(self, child, branch_length=1.0):
        """
        @param child: the child subtree
        @param branch_length: the length of the branch that connects the child subtree
        """
        self.children.append(child)
        child.branch = Branch(self, branch_length)

    def remove_child(self, child):
        child.branch.parent = None
        child.branch = None
        self.children.remove(child)

    def reroot(self):
        """
        Reroot the tree at the current vertex.
        This requires some bookkeeping to keep the correct branch lengths.
        """
        # get the path to the root
        path = []
        v = self
        while v.branch:
            path.append((v.branch.parent, v, v.branch.length))
            v = v.branch.parent
        # reverse parent child relationships on the path
        for parent, child, branch_length in path:
            parent.remove_child(child)
            child.add_child(parent, branch_length)

    def get_unordered_edges(self):
        """
        This function is compatible with the more fundamental unweighted phylogenetic tree.
        @return: a list of ordered vertex pairs
        """
        edges = []
        for v in self.get_vertices():
            if v.branch:
                edges.append((v, v.branch.parent))
                edges.append((v.branch.parent, v))
        return edges

    def get_branches(self):
        return [v.branch for v in self.get_vertices() if v.branch]

    def get_splits(self):
        """
        Get label splits defined by the edges of the tree if self is the root.
        This function calls the base class after some error checking.
        @return: a set of splits
        """
        for branch in self.get_branches():
            if branch.length <= 0:
                raise ValueError('splits on weighted trees are undefined when branch lengths are non-positive')
        return PXVertex.get_splits(self)

    def get_distance_matrix(self):
        """
        This assumes that the labels are integers in [0, N).
        @return: a row major matrix
        """
        labels = self.get_labels()
        n = len(labels)
        if labels != range(n):
            raise ValueError()
        ordered_tips = [None]*n
        for vertex in self.get_labeled_vertices():
            ordered_tips[vertex.label] = vertex
        D = [[0]*n for i in range(n)]
        for i, vertex in enumerate(ordered_tips):
            distance = 0
            stack = [(vertex, vertex.get_neighbor_and_branch_length())]
            while stack:
                source, (sink, distance) = stack.pop()
                if sink.has_label():
                    D[i][sink.label] = distance
                for next_sink, step_distance in sink.get_neighbors_and_branch_lengths():
                    if next_sink is not source:
                        stack.append((sink, (next_sink, distance + step_distance)))
        return D

    def get_weighted_adjacency_matrix(self, id_to_index):
        """
        Weights are reciprocals of branch lengths.
        This function ignores vertex labels.
        @param id_to_index: a map that defines the order of rows in the output matrix
        @return: a row major matrix
        """
        n = len(id_to_index)
        A = [[0]*n for i in range(n)]
        for vertex in self.get_vertices():
            index = id_to_index[id(vertex)]
            for neighbor, blen in vertex.get_neighbors_and_branch_lengths():
                weight = 1.0 / blen
                neighbor_index = id_to_index[id(neighbor)]
                A[index][neighbor_index] = weight
        return A

    def get_newick_substring(self):
        """
        @return: part of a newick string corresponding to the subtree
        """
        if self.has_label():
            s = 'n%d' % self.label
        else:
            s = '(' + ', '.join(child.get_newick_substring() for child in self.children) + ')'
        if self.branch:
            return '%s:%s' % (s, self.branch.length)
        else:
            return s


def set_to_string(my_set):
    """
    @param my_set: a sortable sequence
    @return: a string
    """
    return '{' + ', '.join(str(x) for x in sorted(my_set)) + '}'

def split_to_string(my_split):
    """
    @param my_split: a frozenset of some frozensets
    @return: a string
    """
    return set_to_string([set_to_string(x) for x in my_split])

def list_to_unweighted_tree(term):
    """
    @param term: a label or a list or a list of lists and labels
    @return: the UPXVertex root of an unweighted xtree
    """
    v = UPXVertex()
    if type(term) in (list, tuple):
        for subterm in term:
            v.add_child(list_to_unweighted_tree(subterm))
    else:
        v.label = term
    return v

def list_to_uniformly_weighted_tree(term, branch_length=1.0):
    """
    Each branch will have the same branch length.
    @param term: integer or recursive list
    @param branch_length: the branch length common to all branches
    @return: the WPXVertex root of a weighted xtree
    """
    v = WPXVertex()
    if type(term) in (list, tuple):
        for subterm in term:
            subtree = list_to_uniformly_weighted_tree(subterm, branch_length)
            v.add_child(subtree, branch_length)
    else:
        v.label = term
    return v

def list_to_weighted_tree(term):
    """
    @param term: (integer, weight) pair or recursive list
    @return: the WPXVertex root of a weighted xtree
    """
    v = WPXVertex()
    if type(term) in (list, tuple):
        for subterm, length in term:
            v.add_child(list_to_weighted_tree(subterm), length)
    else:
        v.label = term
    return v

def make_split(a, b):
    """
    @param a: a sequence of hashable values
    @param b: a sequence of hashable values
    @return: a split
    """
    return frozenset((frozenset(a), frozenset(b)))

def splits_to_rf_distance(splits_a, splits_b):
    """
    Compute the Robinson-Foulds distance between two trees given their split sets.
    There is apparently a linear algorithm that calculates this distance,
    but the algorithm used here is worse than linear.
    @param splits_a: nontrivial splits implied by the first xtree
    @param splits_b: nontrivial splits implied by the second xtree
    """
    return 0.5 * len(splits_a ^ splits_b)

def trees_to_rf_distance(root_a, root_b):
    """
    Compute the Robinson-Foulds distance between two trees.
    There is apparently a linear algorithm that calculates this distance,
    but the algorithm used here is worse than linear.
    @param root_a: the root of the first xtree
    @param root_b: the root of the second xtree
    """
    return splits_to_rf_distance(root_a.get_nontrivial_splits(), root_b.get_nontrivial_splits())

class TestXtree(unittest.TestCase):

    def test_nontrivial_splits(self):
        """
        Test the enumeration of nontrivial splits.
        """
        # first test
        root = list_to_unweighted_tree([[0, 1, 2], 3, [4, 5, 6]])
        observed = root.get_nontrivial_splits()
        expected = set((
            make_split((0, 1, 2), (3, 4, 5, 6)),
            make_split((0, 1, 2, 3), (4, 5, 6))))
        self.assertEquals(observed, expected)
        # second test
        root = list_to_unweighted_tree([[0, 1, 2], 3, [[4, 5], [6, 7]]])
        observed = root.get_nontrivial_splits()
        expected = set((
                make_split((0, 1, 2, 3, 4, 5), (6, 7)),
                make_split((0, 1, 2, 3, 6, 7), (4, 5)),
                make_split((0, 1, 2, 3), (4, 5, 6, 7)),
                make_split((0, 1, 2), (3, 4, 5, 6, 7))))
        self.assertEquals(observed, expected)

    def test_distance(self):
        """
        Test the calculation of pairwise distances.
        """
        root = list_to_unweighted_tree([[0, 1, 2], 3, [4, 5, 6]])
        observed = root.get_distance_matrix()
        expected = [
                [0, 2, 2, 3, 4, 4, 4],
                [2, 0, 2, 3, 4, 4, 4],
                [2, 2, 0, 3, 4, 4, 4],
                [3, 3, 3, 0, 3, 3, 3],
                [4, 4, 4, 3, 0, 2, 2],
                [4, 4, 4, 3, 2, 0, 2],
                [4, 4, 4, 3, 2, 2, 0]]
        self.assertEquals(observed, expected)

    def test_newick(self):
        """
        Test the conversion of a tree to a newick string.
        """
        # test the creation of a newick string without branch lengths
        root = list_to_unweighted_tree([[0, 1, 2], 3, [4, 5, 6]])
        observed = root.get_newick_string()
        expected = '((n0, n1, n2), n3, (n4, n5, n6));'
        self.assertEquals(observed, expected)
        # test the creation of a newick string with branch lengths
        root = list_to_uniformly_weighted_tree([[0, 1, 2], 3, [4, 5, 6]], 0.1)
        observed = root.get_newick_string()
        expected = '((n0:0.1, n1:0.1, n2:0.1):0.1, n3:0.1, (n4:0.1, n5:0.1, n6:0.1):0.1);'
        self.assertEquals(observed, expected)

    def test_quartets(self):
        """
        Test the ability to enumerate the quartets implied by a tree.
        """
        root = list_to_unweighted_tree([[0, 1], 2, [3, 4]])
        observed = root.get_quartets()
        expected = set((
            make_split((0, 1), (2, 3)),
            make_split((0, 1), (2, 4)),
            make_split((0, 1), (4, 3)),
            make_split((3, 4), (2, 0)),
            make_split((3, 4), (2, 1))))
        self.assertEquals(observed, expected)

    def test_neighbors_and_branch_lengths(self):
        """
        Test the code that gets the list of neighbors and branch lengths.
        """
        # make the tree
        distance = 2.0
        topo = [0, 1, 2]
        tree = list_to_uniformly_weighted_tree(topo, distance)
        # check the values for the root
        pairs = tree.get_neighbors_and_branch_lengths()
        self.assertEqual(len(pairs), 3)
        for neighbor, blen in pairs:
            self.assertEqual(blen, distance)
        # check the values for the children
        for child in tree.children:
            pairs = child.get_neighbors_and_branch_lengths()
            self.assertEqual(len(pairs), 1)
            for neighbor, blen in pairs:
                self.assertEqual(blen, distance)

    def test_weighted_adjacency(self):
        """
        Test the code that gets the weighted adjacency matrix.
        """
        # define the tree and the id_to_index map
        distance = 2.0
        weight = 0.5
        topo = [0, 1, 2]
        tree = list_to_uniformly_weighted_tree(topo, distance)
        id_to_index = {}
        id_to_index[id(tree.children[0])] = 0
        id_to_index[id(tree.children[1])] = 1
        id_to_index[id(tree.children[2])] = 2
        id_to_index[id(tree)] = 3
        # define the expected weighted adjacency matrix
        expected = (
                (0, 0, 0, weight),
                (0, 0, 0, weight),
                (0, 0, 0, weight),
                (weight, weight, weight, 0))
        # get the observed weighted adjacency matrix
        A = tree.get_weighted_adjacency_matrix(id_to_index)
        observed = tuple(tuple(v) for v in A)
        self.assertEquals(expected, observed)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestXtree)
    unittest.TextTestRunner(verbosity=2).run(suite)

