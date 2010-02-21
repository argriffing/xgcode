"""Check a property of some random graphs.

For each vertex in each random graph,
see if Schur complementing out the vertex results in a graph
whose Fiedler cut is incompatible with the original graph.
Use the Erdos-Renyi probability distribution
over graphs with a given number of vertices.
Edge weights are exponentially distributed.
"""

from StringIO import StringIO
import random

import numpy as np

from SnippetUtil import HandlingError
import Euclid
import MatrixUtil
import BuildTreeTopology
from Form import RadioItem
import Form

def get_form():
    """
    @return: the body of a form
    """
    # define the list of form objects
    form_objects = [
            Form.Integer('ngraphs', 'check this many graphs',
                10, low=1, high=100),
            Form.Integer('nvertices', 'use this many vertices per graph',
                10, low=3, high=20),
            Form.Float('pedge', 'the probability of a vertex pair edge',
                '0.9', low_exclusive=0, high_inclusive=1),
            Form.RadioGroup('cut', 'reduced graph cut criterion', [
                RadioItem('fiedler_cut', 'fiedler cut', True),
                RadioItem('random_cut', 'random cut')])]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object decorated with field values
    @return: a (response_headers, response_text) pair
    """
    out = StringIO()
    # try to make some graphs
    unconnected_count = 0
    invalid_split_count = 0
    valid_split_count = 0
    for graph_index in range(fs.ngraphs):
        G = erdos_renyi(fs.nvertices, fs.pedge)
        if is_connected(G):
            # add interesting edge weights
            add_exponential_weights(G)
            # turn the adjacency matrix into a laplacian matrix
            L = Euclid.adjacency_to_laplacian(G)
            for v in range(fs.nvertices):
                small_index_to_big_index = {}
                for i_small, i_big in enumerate([i for i in range(fs.nvertices) if i != v]):
                    small_index_to_big_index[i_small] = i_big
                # take the schur complement with respect to the given vertex
                L_reduced = get_single_element_schur_complement(L, v)
                assert len(L_reduced) == len(L) - 1
                # get the loadings of the vertices of the reduced graph
                if fs.fiedler_cut:
                    Y_reduced = BuildTreeTopology.laplacian_to_fiedler(L_reduced)
                elif fs.random_cut:
                    Y_reduced = get_random_vector(L_reduced)
                assert len(Y_reduced) == len(L_reduced)
                # expand the fiedler vector with positive and negative valuations for the removed vertex
                found_valid_split = False
                for augmented_loading in (-1.0, 1.0):
                    # get the augmented split vector for this assignment of the removed vertex
                    Y_full = [0]*len(G)
                    for i_reduced, loading in enumerate(Y_reduced):
                        i_big = small_index_to_big_index[i_reduced]
                        Y_full[i_big] = loading
                    Y_full[v] = augmented_loading
                    assert len(Y_full) == len(G)
                    # get the two graphs defined by the split
                    subgraph_a, subgraph_b = list(gen_subgraphs(G, Y_full))
                    # if the subgraphs are both connected then the split is valid
                    if is_connected(subgraph_a) and is_connected(subgraph_b):
                        found_valid_split = True
                # if a valid split was not found then show the matrix
                if found_valid_split:
                    valid_split_count += 1
                else:
                    print >> out, 'Found a matrix that was split incompatibly by a cut of its schur complement!'
                    print >> out, 'matrix:'
                    print >> out, MatrixUtil.m_to_string(G)
                    print >> out, 'index that was removed:', v
                    invalid_split_count += 1
        else:
            unconnected_count += 1
    # show the number of connected and of unconnected graphs
    print >> out, 'this many random graphs were connected:', fs.ngraphs - unconnected_count
    print >> out, 'this many random graphs were not connected:', unconnected_count
    print >> out, 'this many splits were valid:', valid_split_count
    print >> out, 'this many splits were invalid:', invalid_split_count
    # write the result
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue()

def get_single_element_schur_complement(L, index):
    """
    @param L: a laplacian
    @param index: the index of the vertex to be removed
    @return: the laplacian of a graph with one fewer vertex
    """
    # define the number of vertices in the original graph
    n_big = len(L)
    # define the number of vertices in the reduced graph
    n_small = n_big - 1
    # create the reduced graph
    L_reduced = np.zeros((n_small, n_small))
    # define the vertices of the big graph that are also in the small graph
    big_vertices = set(range(n_big)) - set([index])
    for i_small, i_big in enumerate(sorted(big_vertices)):
        for j_small, j_big in enumerate(sorted(big_vertices)):
            L_reduced[i_small][j_small] = L[i_big][j_big]
    # define the side vector to be used for the schur complement
    side = [L[index][i] for i in sorted(big_vertices)]
    corner = L[index][index]
    # subtract a value from each element of the reduced matrix
    for i, side_a in enumerate(side):
        for j, side_b in enumerate(side):
            L_reduced[i][j] -= side_a * side_b / corner
    # return the laplacian of the smaller graph
    return L_reduced

def get_random_vector(L):
    """
    Get a random cut of the graph.
    @param L: a laplacian
    @return: a vector of mostly random positive and negative loadings
    """
    n = len(L)
    # make a vector where each element is negative one or positive one
    Y = [random.choice((-1.0, 1.0)) for i in range(n)]
    # select two elements without replacement and make them negative and positive respectively
    pos_index, neg_index = random.sample(range(n), 2)
    Y[pos_index] = 1.0
    Y[neg_index] = -1.0
    return Y

def gen_subgraphs(G, Y):
    """
    @param G: a symmetric adjacency matrix
    @param Y: a vector of real-valued loadings that defines the split
    @return: the two subgraphs of the graph G defined by the split Y
    """
    assert len(Y) == len(G)
    n = len(Y)
    iplus = [i for i in range(n) if Y[i] > 0]
    iminus = [i for i in range(n) if Y[i] <= 0]
    for vsub in (iplus, iminus):
        nsub = len(vsub)
        gsub = np.zeros((nsub, nsub))
        for small_i, big_i in enumerate(vsub):
            for small_j, big_j in enumerate(vsub):
                gsub[small_i][small_j] = G[big_i][big_j]
        yield gsub

def is_connected(G):
    """
    @param G: a symmetric adjacency matrix
    @return: True if the graph is connected
    """
    n = len(G)
    vertices = set([0])
    observed = set([0])
    while vertices:
        next_vertices = set()
        for v in vertices:
            for i in range(n):
                if i not in observed:
                    if G[v][i]:
                        next_vertices.add(i)
                        observed.add(i)
        vertices = next_vertices
    return observed == set(range(n))

def add_exponential_weights(G):
    """
    Modify the input graph by changing nonzero weights to random exponential weights.
    @param G: a symmetric binary adjacency matrix
    """
    n = len(G)
    for i in range(n):
        for j in range(i):
            if G[i][j]:
                weight = np.random.exponential()
                G[i][j] = weight
                G[j][i] = weight

def erdos_renyi(n, p):
    """
    Sample an Erdos-Renyi random graph.
    @param n: the number of vertices in the graph
    @param p: the probability that two vertices are adjacent
    @return: a symmetric binary adjacency matrix
    """
    G = np.zeros((n, n))
    for i in range(n):
        for j in range(i):
            if random.random() < p:
                G[i][j] = 1.0
                G[j][i] = 1.0
    return G

