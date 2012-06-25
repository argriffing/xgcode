"""
Remove vertices from a graph and change the edge weights.

If edges are like ohms then
the input graph represents a resistor network where
the weight of each edge is the value in ohms of a resistor
on the edge between the two vertices.
The output graph has the same interpretation.
The effective resistance between two vertices in the output graph
should be the same as the effective resistance between the same two
vertices in the input graph.
If edges are like conductances,
then the input and output edges are like reciprocal ohms.
"""

from StringIO import StringIO

import numpy as np

from SnippetUtil import HandlingError
import SnippetUtil
import Util
import MatrixUtil
import iterutils
from Form import RadioItem
import Form
import FormOut
import const

g_edge_data = const.read('20100730u')

def get_form():
    """
    @return: a list of form objects
    """
    # define the edges of the default sparse graph
    edge_lines = Util.get_stripped_lines(g_edge_data.splitlines())
    # define the vertices to be removed by default
    vertices_to_be_removed = ['A', 'a']
    # define the list of form objects
    form_objects = [
            Form.MultiLine('graph', 'sparse graph with edge weights',
                '\n'.join(edge_lines)),
            Form.MultiLine('vertices', 'vertices to be removed',
                '\n'.join(vertices_to_be_removed)),
            Form.RadioGroup('method', 'transformation method', [
                RadioItem('funky', 'some funky method'),
                RadioItem('funky_corrected', 'some funky method (corrected)'),
                RadioItem('ohm', 'edges are ohms'),
                RadioItem('conductance', 'edges are conductances', True)])]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_funky_transformation(edge_triples, name_to_index, reduced_ordered_vertex_names):
    """
    @param edge_triples: TODO
    @param name_to_index: TODO
    @param reduced_ordered_vertex_names: TODO
    @return: reduced edge triples
    """
    # create the graph laplacian
    n = len(name_to_index)
    L = np.zeros((n, n))
    for name_a, name_b, weight in edge_triples:
        a = name_to_index[name_a]
        b = name_to_index[name_b]
        L[a][b] = -weight
        L[b][a] = -weight
        L[a][a] += weight
        L[b][b] += weight
    # get the Moore-Penrose inverse of the laplacian
    L_pinv = np.linalg.pinv(np.array(L))
    # remove rows and columns of specified indices to create a sub matrix
    reduced_name_to_index = dict((name, i) for i, name in enumerate(reduced_ordered_vertex_names))
    reduced_n = len(reduced_ordered_vertex_names)
    reduced_L_pinv = np.zeros((reduced_n, reduced_n))
    for reduced_i in range(reduced_n):
        for reduced_j in range(reduced_n):
            i = name_to_index[reduced_ordered_vertex_names[reduced_i]]
            j = name_to_index[reduced_ordered_vertex_names[reduced_j]]
            reduced_L_pinv[reduced_i][reduced_j] = L_pinv[i][j]
    # get the Moore-Penrose inverse of this reduced matrix
    reduced_L = np.linalg.pinv(reduced_L_pinv)
    # get reduced edge triples
    reduced_edge_triples = []
    epsilon = 0.00000000001
    for i, name_a in enumerate(reduced_ordered_vertex_names):
        for j, name_b in enumerate(reduced_ordered_vertex_names):
            if i < j:
                conductance = -reduced_L[i][j]
                if abs(conductance) > epsilon:
                    triple = (name_a, name_b, conductance)
                    reduced_edge_triples.append(triple)
    return reduced_edge_triples

def get_corrected_funky_transformation(edge_triples, name_to_index, reduced_ordered_vertex_names):
    """
    @param edge_triples: TODO
    @param name_to_index: TODO
    @param reduced_ordered_vertex_names: TODO
    @return: reduced edge triples
    """
    # create the graph laplacian
    n = len(name_to_index)
    L = np.zeros((n, n))
    for name_a, name_b, weight in edge_triples:
        a = name_to_index[name_a]
        b = name_to_index[name_b]
        L[a][b] = -weight
        L[b][a] = -weight
        L[a][a] += weight
        L[b][b] += weight
    # get the Moore-Penrose inverse of the laplacian
    L_pinv = np.linalg.pinv(np.array(L))
    # remove rows and columns of specified indices to create a sub matrix
    reduced_name_to_index = dict((name, i) for i, name in enumerate(reduced_ordered_vertex_names))
    reduced_n = len(reduced_ordered_vertex_names)
    reduced_L_pinv = np.zeros((reduced_n, reduced_n))
    for reduced_i in range(reduced_n):
        for reduced_j in range(reduced_n):
            i = name_to_index[reduced_ordered_vertex_names[reduced_i]]
            j = name_to_index[reduced_ordered_vertex_names[reduced_j]]
            reduced_L_pinv[reduced_i][reduced_j] = L_pinv[i][j]
    # double center the matrix
    reduced_L_pinv = MatrixUtil.double_centered(reduced_L_pinv)
    # get the Moore-Penrose inverse of this reduced matrix
    reduced_L = np.linalg.pinv(reduced_L_pinv)
    # get reduced edge triples
    reduced_edge_triples = []
    epsilon = 0.00000000001
    for i, name_a in enumerate(reduced_ordered_vertex_names):
        for j, name_b in enumerate(reduced_ordered_vertex_names):
            if i < j:
                conductance = -reduced_L[i][j]
                if abs(conductance) > epsilon:
                    triple = (name_a, name_b, conductance)
                    reduced_edge_triples.append(triple)
    return reduced_edge_triples

def get_ohm_transformation(edge_triples, name_to_index, reduced_ordered_vertex_names):
    """
    @param edge_triples: TODO
    @param name_to_index: TODO
    @param reduced_ordered_vertex_names: TODO
    @return: reduced edge triples
    """
    # get the graph laplacian given the resistor network
    n = len(name_to_index)
    L = np.zeros((n, n))
    for name_a, name_b, ohm in edge_triples:
        conductance = 1/ohm
        a = name_to_index[name_a]
        b = name_to_index[name_b]
        L[a][b] = -conductance
        L[b][a] = -conductance
        L[a][a] += conductance
        L[b][b] += conductance
    # get the effective resistance matrix 
    L_pinv = np.linalg.pinv(L) 
    R = np.zeros((n, n)) 
    for i in range(n): 
        for j in range(n): 
            R[i][j] = L_pinv[i][i] + L_pinv[j][j] - L_pinv[i][j] - L_pinv[j][i] 
    # remove rows and columns of specified indices to create a reduced resistance matrix
    reduced_name_to_index = dict((name, i) for i, name in enumerate(reduced_ordered_vertex_names))
    reduced_n = len(reduced_ordered_vertex_names)
    reduced_R = np.zeros((reduced_n, reduced_n))
    for reduced_i in range(reduced_n):
        for reduced_j in range(reduced_n):
            i = name_to_index[reduced_ordered_vertex_names[reduced_i]]
            j = name_to_index[reduced_ordered_vertex_names[reduced_j]]
            reduced_R[reduced_i][reduced_j] = R[i][j]
    # find the laplacian that corresponds to this reduced resistance matrix
    reduced_L = -2*np.linalg.pinv(MatrixUtil.double_centered(reduced_R))
    # get reduced edge triples
    reduced_edge_triples = []
    epsilon = 1e-11
    for i, name_a in enumerate(reduced_ordered_vertex_names):
        for j, name_b in enumerate(reduced_ordered_vertex_names):
            if i < j:
                conductance = -reduced_L[i][j]
                if abs(conductance) > epsilon:
                    ohm = 1/conductance
                    triple = (name_a, name_b, ohm)
                    reduced_edge_triples.append(triple)
    return reduced_edge_triples

def get_conductance_transformation(edge_triples, name_to_index, reduced_ordered_vertex_names):
    """
    @param edge_triples: TODO
    @param name_to_index: TODO
    @param reduced_ordered_vertex_names: TODO
    @return: reduced edge triples
    """
    # get the graph laplacian given the conductance network
    n = len(name_to_index)
    L = np.zeros((n, n))
    for name_a, name_b, conductance in edge_triples:
        a = name_to_index[name_a]
        b = name_to_index[name_b]
        L[a][b] = -conductance
        L[b][a] = -conductance
        L[a][a] += conductance
        L[b][b] += conductance
    # get the effective resistance matrix 
    L_pinv = np.linalg.pinv(L) 
    R = np.zeros((n, n)) 
    for i in range(n): 
        for j in range(n): 
            R[i][j] = L_pinv[i][i] + L_pinv[j][j] - L_pinv[i][j] - L_pinv[j][i] 
    # remove rows and columns of specified indices to create a reduced resistance matrix
    reduced_name_to_index = dict((name, i) for i, name in enumerate(reduced_ordered_vertex_names))
    reduced_n = len(reduced_ordered_vertex_names)
    reduced_R = np.zeros((reduced_n, reduced_n))
    for reduced_i in range(reduced_n):
        for reduced_j in range(reduced_n):
            i = name_to_index[reduced_ordered_vertex_names[reduced_i]]
            j = name_to_index[reduced_ordered_vertex_names[reduced_j]]
            reduced_R[reduced_i][reduced_j] = R[i][j]
    # find the laplacian that corresponds to this reduced resistance matrix
    reduced_L = -2*np.linalg.pinv(MatrixUtil.double_centered(reduced_R))
    # get reduced edge triples
    reduced_edge_triples = []
    epsilon = 1e-11
    for i, name_a in enumerate(reduced_ordered_vertex_names):
        for j, name_b in enumerate(reduced_ordered_vertex_names):
            if i < j:
                conductance = -reduced_L[i][j]
                if abs(conductance) > epsilon:
                    triple = (name_a, name_b, conductance)
                    reduced_edge_triples.append(triple)
    return reduced_edge_triples

def get_response_content(fs):
    # read the edge triples (vertex name, vertex name, edge weight)
    edge_triples = []
    for line in iterutils.stripped_lines(fs.graph.splitlines()):
        string_triple = line.split()
        if len(string_triple) != 3:
            raise HandlingError(
                    'each graph row should have three elements '
                    'but found this line: ' + line)
        triple = string_triple[:2]
        try:
            weight = float(string_triple[2])
        except ValueError as e:
            raise HandlingError(
                    'edge weights should be floating point numbers')
        if weight <= 0:
            raise HandlingError('edge weights should be positive')
        triple.append(weight)
        edge_triples.append(triple)
    # get the set of directed edges to check for redundant or invalid input
    unordered_directed_edges = set()
    for a, b, weight in edge_triples:
        if a == b:
            raise HandlingError(
                    'vertices should not have edges connecting to themselves')
        if (a, b) in unordered_directed_edges:
            raise HandlingError('each edge should be given only once')
        if (b, a) in unordered_directed_edges:
            raise HandlingError(
                    'each edge should be given in only one direction')
        unordered_directed_edges.add((a, b))
    # get the lexicographically ordered list of vertex names
    unordered_vertex_names = set()
    for edge in unordered_directed_edges:
        unordered_vertex_names.update(set(edge))
    ordered_vertex_names = list(sorted(unordered_vertex_names))
    name_to_index = dict(
            (name, i) for i, name in enumerate(ordered_vertex_names))
    n = len(ordered_vertex_names)
    # read the set of vertices that the user wants to remove
    vertex_names_to_remove = set()
    for name in iterutils.stripped_lines(fs.vertices.splitlines()):
        if name in vertex_names_to_remove:
            raise HandlingError(
                    'vertices should be named for removal at most once')
        vertex_names_to_remove.add(name)
    # Assert that the set of vertex names for removal
    # is a subset of the vertex names in the graph.
    weird_names = vertex_names_to_remove - unordered_vertex_names
    if weird_names:
        raise HandlingError(
                'some vertices named for removal '
                'were not found in the graph: ' + str(weird_names))
    # get the ordered list of vertex names that will remain
    reduced_ordered_vertex_names = list(
            sorted(unordered_vertex_names - vertex_names_to_remove))
    # get the laplacian depending on the method
    if fs.funky:
        reduced_edge_triples = get_funky_transformation(
                edge_triples, name_to_index, reduced_ordered_vertex_names)
    elif fs.funky_corrected:
        reduced_edge_triples = get_corrected_funky_transformation(
                edge_triples, name_to_index, reduced_ordered_vertex_names)
    elif fs.ohm:
        reduced_edge_triples = get_ohm_transformation(
                edge_triples, name_to_index, reduced_ordered_vertex_names)
    elif fs.conductance:
        reduced_edge_triples = get_conductance_transformation(
                edge_triples, name_to_index, reduced_ordered_vertex_names)
    # write the reduced edge triples
    out = StringIO()
    for name_a, name_b, weight in reduced_edge_triples:
        print >> out, name_a, name_b, weight
    # write the response
    return out.getvalue()
