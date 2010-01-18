"""Remove vertices from a graph and change the edge weights.

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

import StringIO

from scipy import linalg
import numpy

from SnippetUtil import HandlingError
import SnippetUtil
import Form
import Util
import MatrixUtil

def get_form():
    """
    @return: a list of form objects
    """
    # define the edges of the default sparse graph
    edge_lines = [
            'A a 1',
            'A B 2',
            'A C 2',
            'B D 2',
            'B E 2',
            'C D 2',
            'C E 2',
            'D E 2',
            'a b 3',
            'a c 3',
            'b d 3',
            'b e 3',
            'c d 3',
            'c e 3',
            'd e 3']
    # define the vertices to be removed by default
    vertices_to_be_removed = ['A', 'a']
    # define the list of form objects
    form_objects = [
            Form.MultiLine('graph', 'sparse graph with edge weights', '\n'.join(edge_lines)),
            Form.MultiLine('vertices', 'vertices to be removed', '\n'.join(vertices_to_be_removed)),
            Form.RadioGroup('method', 'transformation method', [
                Form.RadioItem('funky', 'some funky method'),
                Form.RadioItem('funky_corrected', 'some funky method (corrected; edges are like conductances)'),
                Form.RadioItem('ohm', 'edges are like ohms'),
                Form.RadioItem('conductance', 'edges are like conductances', True)])]
    return form_objects

def get_funky_transformation(edge_triples, name_to_index, reduced_ordered_vertex_names):
    """
    @param edge_triples: TODO
    @param name_to_index: TODO
    @param reduced_ordered_vertex_names: TODO
    @return: reduced edge triples
    """
    # create the graph laplacian
    n = len(name_to_index)
    L = numpy.zeros((n, n))
    for name_a, name_b, weight in edge_triples:
        a = name_to_index[name_a]
        b = name_to_index[name_b]
        L[a][b] = -weight
        L[b][a] = -weight
        L[a][a] += weight
        L[b][b] += weight
    # get the Moore-Penrose inverse of the laplacian
    L_pinv = linalg.pinv(numpy.array(L))
    # remove rows and columns of specified indices to create a sub matrix
    reduced_name_to_index = dict((name, i) for i, name in enumerate(reduced_ordered_vertex_names))
    reduced_n = len(reduced_ordered_vertex_names)
    reduced_L_pinv = numpy.zeros((reduced_n, reduced_n))
    for reduced_i in range(reduced_n):
        for reduced_j in range(reduced_n):
            i = name_to_index[reduced_ordered_vertex_names[reduced_i]]
            j = name_to_index[reduced_ordered_vertex_names[reduced_j]]
            reduced_L_pinv[reduced_i][reduced_j] = L_pinv[i][j]
    # get the Moore-Penrose inverse of this reduced matrix
    reduced_L = linalg.pinv(reduced_L_pinv)
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
    L = numpy.zeros((n, n))
    for name_a, name_b, weight in edge_triples:
        a = name_to_index[name_a]
        b = name_to_index[name_b]
        L[a][b] = -weight
        L[b][a] = -weight
        L[a][a] += weight
        L[b][b] += weight
    # get the Moore-Penrose inverse of the laplacian
    L_pinv = linalg.pinv(numpy.array(L))
    # remove rows and columns of specified indices to create a sub matrix
    reduced_name_to_index = dict((name, i) for i, name in enumerate(reduced_ordered_vertex_names))
    reduced_n = len(reduced_ordered_vertex_names)
    reduced_L_pinv = numpy.zeros((reduced_n, reduced_n))
    for reduced_i in range(reduced_n):
        for reduced_j in range(reduced_n):
            i = name_to_index[reduced_ordered_vertex_names[reduced_i]]
            j = name_to_index[reduced_ordered_vertex_names[reduced_j]]
            reduced_L_pinv[reduced_i][reduced_j] = L_pinv[i][j]
    # double center the matrix
    reduced_L_pinv = MatrixUtil.double_centered(reduced_L_pinv)
    # get the Moore-Penrose inverse of this reduced matrix
    reduced_L = linalg.pinv(reduced_L_pinv)
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
    L = numpy.zeros((n, n))
    for name_a, name_b, ohm in edge_triples:
        conductance = 1/ohm
        a = name_to_index[name_a]
        b = name_to_index[name_b]
        L[a][b] = -conductance
        L[b][a] = -conductance
        L[a][a] += conductance
        L[b][b] += conductance
    # get the effective resistance matrix 
    L_pinv = linalg.pinv(L) 
    R = numpy.zeros((n, n)) 
    for i in range(n): 
        for j in range(n): 
            R[i][j] = L_pinv[i][i] + L_pinv[j][j] - L_pinv[i][j] - L_pinv[j][i] 
    # remove rows and columns of specified indices to create a reduced resistance matrix
    reduced_name_to_index = dict((name, i) for i, name in enumerate(reduced_ordered_vertex_names))
    reduced_n = len(reduced_ordered_vertex_names)
    reduced_R = numpy.zeros((reduced_n, reduced_n))
    for reduced_i in range(reduced_n):
        for reduced_j in range(reduced_n):
            i = name_to_index[reduced_ordered_vertex_names[reduced_i]]
            j = name_to_index[reduced_ordered_vertex_names[reduced_j]]
            reduced_R[reduced_i][reduced_j] = R[i][j]
    # find the laplacian that corresponds to this reduced resistance matrix
    reduced_L = -2*linalg.pinv(MatrixUtil.double_centered(reduced_R))
    # get reduced edge triples
    reduced_edge_triples = []
    epsilon = 0.00000000001
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
    L = numpy.zeros((n, n))
    for name_a, name_b, conductance in edge_triples:
        a = name_to_index[name_a]
        b = name_to_index[name_b]
        L[a][b] = -conductance
        L[b][a] = -conductance
        L[a][a] += conductance
        L[b][b] += conductance
    # get the effective resistance matrix 
    L_pinv = linalg.pinv(L) 
    R = numpy.zeros((n, n)) 
    for i in range(n): 
        for j in range(n): 
            R[i][j] = L_pinv[i][i] + L_pinv[j][j] - L_pinv[i][j] - L_pinv[j][i] 
    # remove rows and columns of specified indices to create a reduced resistance matrix
    reduced_name_to_index = dict((name, i) for i, name in enumerate(reduced_ordered_vertex_names))
    reduced_n = len(reduced_ordered_vertex_names)
    reduced_R = numpy.zeros((reduced_n, reduced_n))
    for reduced_i in range(reduced_n):
        for reduced_j in range(reduced_n):
            i = name_to_index[reduced_ordered_vertex_names[reduced_i]]
            j = name_to_index[reduced_ordered_vertex_names[reduced_j]]
            reduced_R[reduced_i][reduced_j] = R[i][j]
    # find the laplacian that corresponds to this reduced resistance matrix
    reduced_L = -2*linalg.pinv(MatrixUtil.double_centered(reduced_R))
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

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # read the edge triples (vertex name, vertex name, edge weight)
    edge_triples = []
    for line in Util.stripped_lines(StringIO.StringIO(fs.graph)):
        string_triple = line.split()
        if len(string_triple) != 3:
            raise HandlingError('each graph edge should have two vertex names and a weight')
        triple = string_triple[:2]
        try:
            weight = float(string_triple[2])
        except ValueError, e:
            raise HandlingError('edge weights should be floating point numbers')
        if weight <= 0:
            raise HandlingError('edge weights should be positive')
        triple.append(weight)
        edge_triples.append(triple)
    # get the set of directed edges to check for redundant or invalid input
    unordered_directed_edges = set()
    for a, b, weight in edge_triples:
        if a == b:
            raise HandlingError('vertices should not have edges connecting to themselves')
        if (a, b) in unordered_directed_edges:
            raise HandlingError('each edge should be given only once')
        if (b, a) in unordered_directed_edges:
            raise HandlingError('each edge should be given in only one direction')
        unordered_directed_edges.add((a, b))
    # get the lexicographically ordered list of vertex names
    unordered_vertex_names = set()
    for edge in unordered_directed_edges:
        unordered_vertex_names.update(set(edge))
    ordered_vertex_names = list(sorted(unordered_vertex_names))
    name_to_index = dict((name, i) for i, name in enumerate(ordered_vertex_names))
    n = len(ordered_vertex_names)
    # read the set of vertices that the user wants to remove
    vertex_names_to_remove = set()
    for name in Util.stripped_lines(StringIO.StringIO(fs.vertices)):
        if name in vertex_names_to_remove:
            raise HandlingError('vertices should be named for removal at most once')
        vertex_names_to_remove.add(name)
    # assert that the set of vertex names for removal is a subset of the vertex names in the graph
    weird_names = vertex_names_to_remove - unordered_vertex_names
    if weird_names:
        raise HandlingError('some vertices named for removal were not found in the graph: ' + str(weird_names))
    # get the ordered list of vertex names that will remain
    reduced_ordered_vertex_names = list(sorted(unordered_vertex_names - vertex_names_to_remove))
    # get the laplacian depending on the method
    if fs.funky:
        reduced_edge_triples = get_funky_transformation(edge_triples, name_to_index, reduced_ordered_vertex_names)
    elif fs.funky_corrected:
        reduced_edge_triples = get_corrected_funky_transformation(edge_triples, name_to_index, reduced_ordered_vertex_names)
    elif fs.ohm:
        reduced_edge_triples = get_ohm_transformation(edge_triples, name_to_index, reduced_ordered_vertex_names)
    elif fs.conductance:
        reduced_edge_triples = get_conductance_transformation(edge_triples, name_to_index, reduced_ordered_vertex_names)
    # write the reduced edge triples
    out = StringIO.StringIO()
    for name_a, name_b, weight in reduced_edge_triples:
        print >> out, name_a, name_b, weight
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()

