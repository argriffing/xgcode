"""Given a newick tree, compare distance and dissimilarity splits.

Distance matrices and dissimilarity matrices are both EDMs
in the mathematical sense.
Eric Stone proved that a vector derived from the distance matrix
will split the tree,
and here I am trying to show that a vector similarly
derived from a dissimilarity matrix does not necessarily
split the tree in a way that defines a unique branch.
"""

import StringIO
import math

import numpy
from scipy import linalg

from SnippetUtil import HandlingError
import Form
import MatrixUtil
import NewickIO
import FelTree

def distance_to_dissimilarity(d, k):
    """
    @param d: the additive distance between points on a tree
    @param k: the residue alphabet size
    @return: the expected proportion of sites that do not match
    """
    b = float(k) / float(k-1)
    return (1 - math.exp(-b*d)) / b

def get_principal_eigenvector(M):
    """
    @param M: a 2d numpy array representing a matrix
    @return: the principal eigenvector of M
    """
    eigenvalues, eigenvector_transposes = linalg.eigh(M)
    eigenvectors = eigenvector_transposes.T
    eigensystem = [(abs(w), w, v.tolist()) for w, v in zip(eigenvalues, eigenvectors)]
    sorted_eigensystem = list(reversed(sorted(eigensystem)))
    sorted_abs_eigenvalues, sorted_eigenvalues, sorted_eigenvectors = zip(*sorted_eigensystem)
    principal_eigenvector = sorted_eigenvectors[0]
    return principal_eigenvector

def get_principal_coordinate(D):
    """
    Return the principal eigenvector of the pseudoinverse of the laplacian.
    @param D: a distance matrix
    @return: the principal coordinate
    """
    L_pinv = -0.5 * MatrixUtil.double_centered(D)
    return get_principal_eigenvector(L_pinv)

def get_form():
    """
    @return: a list of form objects
    """
    # this counterexample uses symmetric and long branch lengths
    unused_counterexample = '(a:2, b:5, (c:5, d:2):1);'
    # this counterexample uses asymmetric and shorter branch lengths
    counterexample = '(a:.5, b:2.5, (c:2.7, d:.6):.3);'
    default_tree_string = counterexample
    tree = NewickIO.parse(default_tree_string, FelTree.NewickTree)
    formatted_tree_string = NewickIO.get_narrow_newick_string(tree, 60)
    # return the form objects
    return [Form.MultiLine('tree', 'newick tree with branch lengths', formatted_tree_string)]

def get_response(fs):
    """
    @param fs: a FieldStorage object decorated with field values
    @return: a (response_headers, response_text) pair
    """
    # arbitrarily define the size of the alphabet
    k = 4
    # define the response
    out = StringIO.StringIO()
    # get the tree
    tree = NewickIO.parse(fs.tree, FelTree.NewickTree)
    # define the order of the tip names
    ordered_tip_names = list(sorted(node.get_name() for node in tree.gen_tips()))
    n = len(ordered_tip_names)
    # get the matrix of pairwise distances among the tips
    D = numpy.array(tree.get_distance_matrix(ordered_tip_names))
    D_vector = get_principal_coordinate(D)
    # get the dissimilarity matrix from the distance matrix
    dissimilarity = numpy.array([[distance_to_dissimilarity(d, k) for d in row] for row in D])
    dissimilarity_vector = get_principal_coordinate(dissimilarity)
    # get the principal coordinates of the distance-like matrices
    print >> out, 'original distance matrix:'
    print >> out, MatrixUtil.m_to_string(D)
    print >> out
    print >> out, 'projections onto the principal coordinate using the original distance matrix:'
    for name, value in zip(ordered_tip_names, D_vector):
        print >> out, '\t'.join((name, str(value)))
    print >> out
    print >> out, 'dissimilarity matrix:'
    print >> out, MatrixUtil.m_to_string(dissimilarity)
    print >> out
    print >> out, 'projections onto the principal coordinate using the dissimilarity matrix:'
    for name, value in zip(ordered_tip_names, dissimilarity_vector):
        print >> out, '\t'.join((name, str(value)))
    print >> out
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
