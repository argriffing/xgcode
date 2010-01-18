"""Consider the graph laplacian matrix for a phylogenetic tree.
"""

import StringIO
import math

from scipy import linalg
import scipy
import numpy

from SnippetUtil import HandlingError
import Newick
import Form

def gen_ancestors(node):
    """
    @param node: a node in a newick tree
    @return: the ordered list of ancestors from the root to the node
    """
    ancestors = []
    current = node
    while current:
        ancestors.append(current)
        current = current.parent
    return reversed(ancestors)

def get_mrca(node_a, node_b):
    """
    Get the most recent common ancestral node on the tree.
    @param node_a: a node in a newick tree
    @param node_b: another node in a newick tree
    @return: the ancestor of a and b that is furthest from the root
    """
    mrca = None
    for pa, pb in zip(gen_ancestors(node_a), gen_ancestors(node_b)):
        if pa is pb:
            mrca = pa
        else:
            break
    return mrca

def get_distance(node_a, node_b):
    """
    @param node_a: a node in a newick tree
    @param node_b: another node in a newick tree
    @return: the branch length between the nodes
    """
    # get the most recent common ancestral node of a and b
    mrca = get_mrca(node_a, node_b)
    # initialize the branch length between the two nodes
    blen = 0
    # add the length to node a
    current = node_a
    while current is not mrca:
        blen += current.blen
        current = current.parent
    # add the length to node b
    current = node_b
    while current is not mrca:
        blen += current.blen
        current = current.parent
    # return the branch length between the nodes
    return blen

def get_form():
    """
    @return: the body of a form
    """
    default_tree_string = '(a:1, (b:2, d:5):1, c:4);'
    return [Form.MultiLine('tree', 'newick tree', default_tree_string)]

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # get the tree
    tree = Newick.parse(fs.tree, Newick.NewickTree)
    tree.assert_valid()
    # get the lexicographically ordered tip names
    states = list(sorted(node.name for node in tree.gen_tips()))
    # start to prepare the reponse
    out = StringIO.StringIO()
    # create the dictionary distance matrix
    dictionary_distance_matrix = {}
    for ta in tree.gen_tips():
        for tb in tree.gen_tips():
            key = (ta.name, tb.name)
            distance = get_distance(ta, tb)
            dictionary_distance_matrix[key] = distance
    # show the distance matrix
    print >> out, 'distance matrix:'
    for a in states:
        for b in states:
            print >> out, a, b, dictionary_distance_matrix[(a, b)]
    print >> out, ''
    # create the off diagonals of a rate matrix from the distance matrix
    unnormalized_rate_matrix = {}
    for key, distance in dictionary_distance_matrix.items():
        a, b = key
        if a == b:
            rate = 0
        else:
            rate = 1.0 / distance
        unnormalized_rate_matrix[key] = rate
    # create the diagonals of the rate matrix
    for a in states:
        row_rate = sum(unnormalized_rate_matrix[(a, b)] for b in states)
        unnormalized_rate_matrix[(a, a)] = -row_rate
    # show the unnormalized rate matrix
    print >> out, 'unnormalized rate matrix:'
    for a in states:
        for b in states:
            print >> out, a, b, unnormalized_rate_matrix[(a, b)]
    print >> out, ''
    """
    # normalize the off diagonals and put negative ones on the diagonal
    rate_matrix = dict(unnormalized_rate_matrix.items())
    for a in states:
        normalizing_factor = -1.0 / rate_matrix[(a, a)]
        for b in states:
            rate_matrix[(a, b)] *= normalizing_factor
    # show the rate matrix
    print >> out, 'normalized rate matrix:'
    for a in states:
        for b in states:
            print >> out, a, b, rate_matrix[(a, b)]
    print >> out, ''
    """
    # create the numpy rate matrix
    row_major = []
    for a in states:
        row = [unnormalized_rate_matrix[(a, b)] for b in states]
        row_major.append(row)
    L = numpy.array(row_major)
    D = scipy.diag([math.sqrt(-1/unnormalized_rate_matrix[(a, a)]) for a in states])
    numpy_matrix = numpy.dot(D, numpy.dot(L, D))
    print >> out, 'symmetrically normalized matrix:'
    print >> out, numpy_matrix
    print >> out, ''
    # show the eigendecomposition of the matrix
    w, vl, vr = linalg.eig(numpy_matrix, left=True, right=True)
    print >> out, 'eigenvalues:'
    print >> out, w
    print >> out, 'left eigenvectors:'
    print >> out, vl
    print >> out, 'right eigenvectors:'
    print >> out, vr
    # get the eigenvalues sorted by absolute value
    ordered_eigenvalue_info = list(sorted((abs(x), i) for i, x in enumerate(w)))
    # get the index of the eigenvector with the second smallest absolute value
    fiedler_eigenvalue, fiedler_eigenvalue_index = ordered_eigenvalue_info[1]
    fiedler_vector = vl.T[fiedler_eigenvalue_index]
    print >> out, 'second smallest absolute value of any eigenvalue:'
    print >> out, fiedler_eigenvalue
    print >> out, 'corresponding fiedler vector:'
    print >> out, fiedler_vector
    print >> out, 'corresponding partition:'
    print >> out, [name for name, value in zip(states, fiedler_vector) if value < 0]
    print >> out, [name for name, value in zip(states, fiedler_vector) if value >= 0]
    # get the laplacian according to the method of balaji and bapat
    # create the numpy distance matrix
    row_major = []
    for a in states:
        row = [dictionary_distance_matrix[(a, b)] for b in states]
        row_major.append(row)
    numpy_distance_matrix = numpy.array(row_major)
    # get the inverse of the numpy distance matrix,
    # using the notation of balaji and bapat
    E_inv = linalg.inv(numpy_distance_matrix)
    ones = numpy.ones((len(states), len(states)))
    numerator = numpy.dot(E_inv, numpy.dot(ones, E_inv))
    denominator = sum(sum(E_inv))
    R = E_inv - numerator / denominator
    print >> out, 'R:'
    print >> out, R
    print >> out, 'sum(R):'
    print >> out, sum(R)
    print >> out, 'eigenvalues of R:'
    w, v = linalg.eig(R)
    print >> out, w
    # get the eigenvalues sorted by absolute value
    ordered_eigenvalue_info = list(sorted((abs(x), i) for i, x in enumerate(w)))
    # get the index of the eigenvector with the second smallest absolute value
    fiedler_eigenvalue, fiedler_eigenvalue_index = ordered_eigenvalue_info[1]
    fiedler_vector = v.T[fiedler_eigenvalue_index]
    print >> out, 'second smallest absolute value of any eigenvalue:'
    print >> out, fiedler_eigenvalue
    print >> out, 'corresponding fiedler vector:'
    print >> out, fiedler_vector
    print >> out, 'corresponding partition:'
    print >> out, [name for name, value in zip(states, fiedler_vector) if value < 0]
    print >> out, [name for name, value in zip(states, fiedler_vector) if value >= 0]
    # try an alternative formulation for R
    n = len(states)
    P = numpy.eye(n) - numpy.ones((n, n)) / n
    R_inv = numpy.dot(P, numpy.dot(numpy_distance_matrix, P))
    print >> out, 'alternative formulation of R inverse:'
    print >> out, R_inv
    w, v = linalg.eigh(R_inv)
    ordered_eigenvalue_info = list(sorted((abs(x), i) for i, x in enumerate(w)))
    fiedler_eigenvalue, fiedler_eigenvalue_index = ordered_eigenvalue_info[1]
    fiedler_vector = v.T[fiedler_eigenvalue_index]
    print >> out, 'eigensystem:'
    print >> out, w
    print >> out, v
    print >> out, 'corresponding fiedler vector:'
    print >> out, fiedler_vector
    print >> out, 'alternative formulation of R:'
    R = linalg.pinv(R_inv)
    print >> out, R
    w, v = linalg.eigh(R)
    ordered_eigenvalue_info = list(sorted((abs(x), i) for i, x in enumerate(w)))
    fiedler_eigenvalue, fiedler_eigenvalue_index = ordered_eigenvalue_info[1]
    fiedler_vector = v.T[fiedler_eigenvalue_index]
    print >> out, 'eigensystem:'
    print >> out, w
    print >> out, v
    print >> out, 'corresponding fiedler vector:'
    print >> out, fiedler_vector
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
