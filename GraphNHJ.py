"""
This module uses graph representations for neighborhood joining algorithms.
The name GraphNHJ is a shortened name for graph based neighborhood joining.
"""

from StringIO import StringIO
import unittest

import numpy
from scipy import linalg
import pygraphviz
#import networkx
#import matplotlib.pyplot

import Newick
import NewickIO
import FelTree
import NeighborJoining
import MatrixUtil
import Clustering


def get_augmented_gower_selection(D):
    """
    Do a spectral sign split with neighbor joining fallback.
    The first choice is to return indices corresponding to
    positive elements of the dominant eigenvector of the gower matrix.
    If this defines a degenerate bipartition,
    then neighbor joining is used as a fallback.
    @param D: a distance matrix
    @return: the set of selected indices
    """
    n = len(D)
    if n < 4:
        raise ValueError('expected a distance matrix with at least four rows')
    # get the gower matrix
    G = MatrixUtil.double_centered(numpy.array(D))
    # get the dominant eigenvector
    eigenvalues, eigenvector_transposes = linalg.eigh(G)
    eigenvectors = eigenvector_transposes.T
    dominant_value, dominant_vector = max((abs(w), v) for w, v in zip(eigenvalues, eigenvectors))
    # get the bipartition defined by the dominant eigenvector
    selection = set(i for i, x in enumerate(dominant_vector) if x > 0)
    complement = set(range(n)) - selection
    # if the bipartition is degenerate then resort to neighbor joining
    if min(len(selection), len(complement)) < 2:
        selection = set(NeighborJoining.get_neighbors(D))
    return selection


class Tree:
    """
    This is an unrooted tree.
    Its representation is really more like that of a graph though.
    """

    def __init__(self, D):
        # create the root node that contains the distance matrix
        self.root = Node()
        # initialize the root node
        n = len(D)
        self.root.D = D
        self.root.index_map = dict(zip(range(n), range(n)))
        # initialize a source of unique index values for internal nodes
        self.next_index = n

    def _use_next_index(self):
        """
        Calling this function has a side effect of incrementing the next index.
        @return: the next number to be used as an internal vertex identifier
        """
        index = self.next_index
        self.next_index += 1
        return index

    def get_decomposible_node(self):
        """
        Find a node whose distance matrix has at least four rows.
        @return: None or a decomposible node
        """
        for node in self.gen_nodes():
            if node.get_dimension() >= 4:
                return node
        return None

    def gen_nodes(self):
        """
        Yield each node in the tree.
        """
        links = []
        if self.root:
            links.append(Link(None, self.root, None))
        while links:
            next_links = []
            for link in links:
                yield link.get_sink()
                next_links.extend(list(link.gen_next_links()))
            links = next_links

    def decompose(self, compound_node, splitter):
        """
        Break a node into two nodes, preserving connectivity.
        If the original distance matrix satisfies a tree metric,
        then each decomposition step preserves geodesic distances between taxa.
        @param compound_node: the node to be decomposed
        @param splitter: a function that gets a bipartition given a distance matrix
        """
        D = compound_node.get_distance_matrix()
        n = compound_node.get_dimension()
        # the selection and complement are sets of local indices
        selection = splitter(D)
        complement = set(range(n)) - selection
        ordered_selection = list(sorted(selection))
        ordered_complement = list(sorted(complement))
        # get the vector of distances to a new vertex
        v = get_split_branch_midpoint_distances(D, selection, complement)
        internal_index = self._use_next_index()
        # create the two new nodes
        node_a = Node()
        node_b = Node()
        new_nodes = [node_a, node_b]
        for node, indices in zip(new_nodes, (ordered_selection, ordered_complement)):
            # add inherited indices to the new index map
            for i, i_local in enumerate(indices):
                node.index_map[i] = compound_node.index_map[i_local]
            # add the index of the new vertex to the new index map
            node.index_map[len(indices)] = internal_index
            # create the distance matrix for the new node
            for i, i_local in enumerate(indices):
                row = []
                for j_local in indices:
                    row.append(D[i_local][j_local])
                row.append(v[i_local])
                node.D.append(row)
            row = []
            for j_local in indices:
                row.append(v[j_local])
            row.append(0)
            node.D.append(row)
            # identify the links for which the new node should be the source
            relevant_out_links = [link for link in compound_node.links if link.index in node.index_map.values()]
            # identify the links for which the new node should be the sink
            relevant_in_links = [link.get_back_link() for link in relevant_out_links]
            # refer to the new node instead of the old node in all relevant links
            for link in relevant_out_links:
                link.source = node
            for link in relevant_in_links:
                link.sink = node
            # define the set of links from the new node
            node.links = relevant_out_links
        # link the two new nodes together
        node_a.links.append(Link(node_a, node_b, internal_index))
        node_b.links.append(Link(node_b, node_a, internal_index))
        # if the compound node is the root then change the root to a new node
        if self.root is compound_node:
            self.root = node_a
        # assert that the old node is not reachable
        for node in self.gen_nodes():
            if node is compound_node:
                raise RuntimeError('failed to remove the decomposed node from the graph')
        # assert that the source of each link from a node is correct
        for node in self.gen_nodes():
            for link in node.links:
                if link.source is not node:
                    raise RuntimeError('found a link whose source was invalid')

    def get_pygraphviz_graph(self, node_labels):
        """
        Return a pygraphviz representation of the tree and its fully connected internal components.
        @param node_labels: labels of the nodes in the same order as in the original distance matrix
        @return: a pygraphviz graph object
        """
        # create the initial graph
        G = pygraphviz.AGraph()
        # add labeled edges
        for node in self.gen_nodes():
            indices = set(node.index_map.values())
            for local_i, global_i in node.index_map.items():
                for local_j, global_j in node.index_map.items():
                    if local_i < local_j:
                        # add a labeled edge
                        edge_label = str(node.D[local_i][local_j])
                        G.add_edge(global_i, global_j, label=edge_label)
        # label all nodes with a default label
        for pygraphviz_node in G.nodes():
            pygraphviz_node.attr['label'] = ' '
        # informatively label nodes that represent taxa
        for i, node_label in enumerate(node_labels):
            pygraphviz_node = G.get_node(i)
            pygraphviz_node.attr['label'] = node_label
        # return the graph
        return G

#    def get_networkx_graph(self):
#        """
#        Return a networkx representation of the tree and its fully connected internal components.
#        @param ordered_labels: labels of the nodes in the same order as in the original distance matrix
#        @return: a networkx graph object
#        """
#        G = networkx.Graph()
#        G.add_nodes_from(range(self.next_index))
#        for node in self.gen_nodes():
#            indices = set(node.index_map.values())
#            for i in indices:
#                for j in indices:
#                    if i < j:
#                        G.add_edge(i, j)
#        return G


class Link:
    """
    This is a link between nodes of the unrooted tree.
    """

    def __init__(self, source, sink, index):
        """
        @param source: the source node object
        @param sink: the target node object
        @param index: the index of the vertex shared by the two nodes
        """
        self.source = source
        self.sink = sink
        self.index = index

    def get_index(self):
        return self.index

    def get_sink(self):
        return self.sink

    def get_source(self):
        return self.source

    def get_back_link(self):
        for link in self.sink.links:
            if link.sink is self.source:
                return link

    def gen_next_links(self):
        for link in self.sink.links:
            if link.sink is not self.source:
                yield link


class Node:
    """
    This is a node in an unrooted tree.
    """

    def __init__(self):
        # links to other nodes
        self.links = []
        # the distance matrix for the node
        self.D = []
        # a map from local to global distance matrix indices
        self.index_map = {}

    def get_distance_matrix(self):
        """
        @return: the local distance matrix
        """
        return self.D

    def get_dimension(self):
        """
        @return: the number of rows in the local distance matrix
        """
        return len(self.get_distance_matrix())


def get_split_branch_arbitrary_distances(D, selection, complement):
    """
    Get distances to a translation of a point on the branch defined by the split.
    @param D: a distance matrix
    @param selection: one of the clusters, represented by a set of indices
    @param complement: the other cluser
    @return: the distance from each vertex to a point related to the branch defined by the split
    """
    # get the number of vertices
    n = len(selection | complement)
    # get the expected distance to a point in the other cluster
    E = [0]*n
    for A, B in ((selection, complement), (complement, selection)):
        for i in A:
            E[i] = sum(D[i][j] for j in B)/len(B)
    # get the mean of E
    E_bar = sum(E)/n
    # get the vector of distances to a virtual point that may or may not be on the tree
    v = [0]*n
    for A, B in ((selection, complement), (complement, selection)):
        for i in A:
            v[i] = E[i] - len(A)*E_bar/n
    return v

def get_split_branch_length(D, selection, complement):
    """
    Get the length of the branch defined by the bipartition.
    @param D: a distance matrix
    @param selection: one of the clusters, represented by a set of indices
    @param complement: the other cluser
    @return: the distance from each vertex to the root of the other subtree
    """
    # get the vector of distances to a virtual point that may or may not be on the tree
    v = get_split_branch_arbitrary_distances(D, selection, complement)
    # get two quantities that should sum to twice the length of the branch defined by the split
    a = min(v[i] + v[j] - D[i][j] for i in selection for j in selection)
    b = min(v[i] + v[j] - D[i][j] for i in complement for j in complement)
    branch_length = (a + b) / 2
    # return the branch length
    return branch_length

def get_split_branch_midpoint_distances(D, selection, complement):
    """
    Get distances to the midpoint of the the branch defined by the split.
    @param D: a distance matrix
    @param selection: one of the clusters, represented by a set of indices
    @param complement: the other cluser
    @return: the distance from each vertex to the midpoint of the branch defined by the split
    """
    n = len(selection | complement)
    # get the vector of distances to a virtual point that may or may not be on the tree
    v = get_split_branch_arbitrary_distances(D, selection, complement)
    # get two quantities that should sum to twice the length of the branch defined by the split
    a = min(v[i] + v[j] - D[i][j] for i in selection for j in selection)
    b = min(v[i] + v[j] - D[i][j] for i in complement for j in complement)
    # get the distance to the midpoint of the split
    nv = [0]*n
    for i in selection:
        nv[i] = v[i] + (b - a)/4
    for i in complement:
        nv[i] = v[i] + (a - b)/4
    return nv


class TestGraphNHJ(unittest.TestCase):

    def test_placeholder(self):
        #FIXME
        pass


def get_sample_treegraph():
    """
    This is used for testing graph drawing.
    I would like to test the drawing without actually creating an image file.
    @return: a Tree object and its ordered labels
    """
    tree_string = '((a:1, b:2):1, (c:1, d:2):3, (e:2, f:3):5);'
    feltree = NewickIO.parse(tree_string, FelTree.NewickTree)
    ordered_labels = list(node.get_name() for node in feltree.gen_tips())
    D = feltree.get_distance_matrix(ordered_labels)
    return Tree(D), ordered_labels

def example_draw_pygraphviz():
    """
    This example shows how to draw the graph using pygraphviz.
    """
    treegraph, ordered_labels = get_sample_treegraph()
    for i in range(10):
        print 'iteration', i
        # get the networkx representation of the graph
        G = treegraph.get_pygraphviz_graph(ordered_labels)
        # FIXME show the number of nodes
        print len(list(G.nodes())), 'nodes'
        # try to prevent the nodes from overlapping when they are drawn
        G.graph_attr['overlap'] = 'false'
        # draw the figure
        G.layout()
        G.draw('graph-%d.png' % i)
        # decompose another node if possible
        node = treegraph.get_decomposible_node()
        if node:
            treegraph.decompose(node, get_augmented_gower_selection)
        else:
            break

def example_draw_networkx():
    """
    This example shows how to draw the graph using networkx.
    """
    treegraph, ordered_labels = get_sample_treegraph()
    for i in range(10):
        print 'iteration', i
        # get the networkx representation of the graph
        G = graphtree.get_networkx_graph()
        # clear the figure
        matplotlib.pyplot.clf()
        # draw the figure
        networkx.draw(G)
        # save the figure
        matplotlib.pyplot.savefig('graph-%d.png' % i)
        # decompose another node if possible
        node = treegraph.get_decomposible_node()
        if node:
            treegraph.decompose(node, get_augmented_gower_selection)
        else:
            break

if __name__ == '__main__':
    """
    suite = unittest.TestLoader().loadTestsFromTestCase(TestGraphNHJ)
    unittest.TextTestRunner(verbosity=2).run(suite)
    """
    example_draw_pygraphviz()

