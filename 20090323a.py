"""Given a newick tree, investigate the spherical euclidean distance matrix.

Given a newick tree,
investigate the corresponding spherical euclidean distance matrix.
"""

from StringIO import StringIO
import math

import numpy as np

from SnippetUtil import HandlingError
import Form
import FormOut
import MatrixUtil
import NewickIO
import FelTree
import const

g_default_string = const.read('20100730q')

def get_form():
    """
    @return: a list of form objects
    """
    # define the default tree string
    tree = NewickIO.parse(g_default_string, FelTree.NewickTree)
    formatted_tree_string = NewickIO.get_narrow_newick_string(tree, 60)
    # return the form objects
    form_objects = [
            Form.MultiLine('tree', 'newick tree with branch lengths',
                formatted_tree_string)]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_matrix_square_root(M, epsilon=1e-14):
    """
    Note that the square root is assumed to have rank one less than that of M.
    Eigenvalues less than epsilon are set to zero.
    @param M: a nxn numpy array
    @return: a nxn numpy array X such that np.dot(X, X.T) == M
    """
    eigenvalues, eigenvector_transposes = np.linalg.eigh(M)
    eigenvalues = [(v if abs(v) > epsilon else 0) for v in eigenvalues]
    eigenvectors = eigenvector_transposes.T
    eigensystem = [(abs(w), w, v.tolist()) for w, v in zip(eigenvalues, eigenvectors)]
    sorted_eigensystem = list(reversed(sorted(eigensystem)))
    sorted_abs_eigenvalues, sorted_eigenvalues, sorted_eigenvectors = zip(*sorted_eigensystem[:-1])
    # for each coordinate, get the value corresponding to each point
    coordinate_array = []
    for value, vector in zip(sorted_eigenvalues, sorted_eigenvectors):
        if value < 0:
            raise HandlingError('negative eigenvalue: %f' % value)
        row = [math.sqrt(value) * v for v in vector]
        coordinate_array.append(row)
    X = np.array(zip(*coordinate_array))
    return X

def euclidean_distance(a, b):
    """
    @param a: a list of coordinates defining a point in euclidean space
    @param b: a list of coordinates defining a point in euclidean space
    """
    distance_squared = sum((x-y)**2 for x, y in zip(a, b))
    return math.sqrt(distance_squared)

def magnitude(v):
    """
    @param v: a list of coordinates defining a point in euclidean space
    """
    distance_squared = sum(x**2 for x in v)
    return math.sqrt(distance_squared)

def column_vector_to_list(v):
    return [row[0] for row in v]

def list_to_column_vector(arr):
    return np.array([[element] for element in arr])

def get_response(fs):
    """
    @param fs: a FieldStorage object decorated with field values
    @return: a (response_headers, response_text) pair
    """
    # define the response
    out = StringIO()
    # get the tree
    tree = NewickIO.parse(fs.tree, FelTree.NewickTree)
    # define the order of the tip names
    ordered_tip_names = list(sorted(node.get_name() for node in tree.gen_tips()))
    n = len(ordered_tip_names)
    # define a conformant column vector of ones and a conformant identity matrix
    e = np.array([[1.0] for i in range(n)])
    I = np.identity(n)
    # get the matrix of pairwise distances among the tips
    D = np.array(tree.get_distance_matrix(ordered_tip_names))
    print >> out, 'original distance matrix:'
    print >> out, MatrixUtil.m_to_string(D)
    print >> out
    print >> out, 'determinant of the original distance matrix:'
    print >> out, np.linalg.det(D)
    print >> out
    D_inv = np.linalg.inv(D)
    # define an s column vector whose elements are positive weights that sum to one
    D_inv_sum = sum(sum(row) for row in D_inv)
    print >> out, 'sum of elements of the inverse distance matrix:'
    print >> out, D_inv_sum
    print >> out
    print >> out, 'square root of one half the inverse of the sum of elements of the inverse distance matrix:'
    print >> out, math.sqrt(.5 * (1/D_inv_sum))
    print >> out
    s = np.dot(D_inv, e) / D_inv_sum
    print >> out, 'the weight vector:'
    print >> out, s
    print >> out
    print >> out, 'the normalized vector of squared weights:'
    squared_elements = [x*x for x in column_vector_to_list(s)]
    print >> out, list_to_column_vector([x/sum(squared_elements) for x in squared_elements])
    print >> out
    print >> out, 'the normalized vector of square root weights:'
    sqrt_elements = [math.sqrt(x) for x in column_vector_to_list(s)]
    print >> out, list_to_column_vector([x/sum(sqrt_elements) for x in sqrt_elements])
    print >> out
    # define the F matrix
    F_left = I - np.dot(e, s.T)
    F_right = I - np.dot(s, e.T)
    F = -0.5 * np.dot(F_left, np.dot(D, F_right))
    # take the matrix square root of the F matrix
    X = get_matrix_square_root(F)
    print >> out, 'euclidean points where the origin is the circumcenter:'
    for row in X:
        print >> out, row
    print >> out
    print >> out, 'distances from the circumcenter (should all be the same):'
    for row in X:
        print >> out, magnitude(row)
    print >> out
    # get the implied distance matrix
    D_implied = np.zeros((n,n))
    for i, row_a in enumerate(X):
        for j, row_b in enumerate(X):
            D_implied[i][j] = euclidean_distance(row_a, row_b)**2
    print >> out, 'implied distance matrix:'
    print >> out, MatrixUtil.m_to_string(D_implied)
    print >> out
    # now set the origin to the center of mass
    F_left = I - np.dot(e, e.T) / np.dot(e.T, e)
    F_right = I - np.dot(e, e.T) / np.dot(e.T, e)
    F = -0.5 * np.dot(F_left, np.dot(D, F_right))
    # take the matrix square root of the F matrix
    X = get_matrix_square_root(F)
    print >> out, 'euclidean points where the origin is the center of mass:'
    for row in X:
        print >> out, row
    print >> out
    # show the center of mass
    print >> out, 'center of mass (should be the origin):'
    center = [sum(values)/len(values) for values in X.T]
    print >> out, center
    print >> out
    # get the implied distance matrix
    D_implied = np.zeros((n,n))
    for i, row_a in enumerate(X):
        for j, row_b in enumerate(X):
            D_implied[i][j] = euclidean_distance(row_a, row_b)**2
    print >> out, 'implied distance matrix:'
    print >> out, MatrixUtil.m_to_string(D_implied)
    print >> out
    print >> out, 'normalized square roots of diagonal of -(1/2)HDH:'
    sqrts_of_diagonal = [math.sqrt(F[i][i]) for i in range(n)]
    normalized_sqrts_of_diagonal = [d / sum(sqrts_of_diagonal) for d in sqrts_of_diagonal]
    for d in normalized_sqrts_of_diagonal:
        print >> out, d
    print >> out
    print >> out, 'distances from the center of mass:'
    distances_from_origin = [magnitude(row) for row in X]
    for distance in distances_from_origin:
        print >> out, distance
    print >> out
    print >> out, 'normalized distances from the center of mass:'
    normalized_distances = [distance / sum(distances_from_origin) for distance in distances_from_origin]
    for distance in normalized_distances:
        print >> out, distance
    print >> out
    print >> out, 'squared distances from the center of mass:'
    squared_distances_from_origin = [magnitude(row)**2 for row in X]
    for distance in squared_distances_from_origin:
        print >> out, distance
    print >> out
    print >> out, 'normalized squared distances from the center of mass:'
    normalized_distances = [distance / sum(squared_distances_from_origin) for distance in squared_distances_from_origin]
    for distance in normalized_distances:
        print >> out, distance
    print >> out
    print >> out, 'L1 distances from the center of mass:'
    L1_distances_from_origin = [sum(abs(x) for x in row) for row in X]
    for distance in L1_distances_from_origin:
        print >> out, distance
    print >> out
    print >> out, 'normalized L1 distances from the center of mass:'
    normalized_L1_distances = [distance / sum(L1_distances_from_origin) for distance in L1_distances_from_origin]
    for distance in normalized_L1_distances:
        print >> out, distance
    print >> out
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
