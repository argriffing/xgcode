"""Given an NxN distance matrix and a split, get two smaller distance matrices.

Given D NxN and a split,
get a (M+1)x(M+1) and a (N-M+1)x(N-M+1) distance matrix.
Given a distance matrix relating tips of the tree,
a conformant list of ordered tip names,
and a list of selected tip names,
use a neighbor-joining-like method to compute
the two new distance matrices,
including nodes that represent the other side of the tree.
"""

from StringIO import StringIO
import itertools

import numpy as np

from SnippetUtil import HandlingError
import Util
import MatrixUtil
import Form
import SchurAlgebra
import BuildTreeTopology

def get_form():
    """
    @return: the body of a form
    """
    # Define the default distance matrix, the ordered labels,
    # and the selected labels.
    D = np.array([
        [ 0,  8, 15, 11, 13,  5],
        [ 8,  0, 11, 15, 17,  9],
        [15, 11,  0, 22, 24, 16],
        [11, 15, 22,  0,  4, 10],
        [13, 17, 24,  4,  0, 12],
        [ 5,  9, 16, 10, 12,  0]], dtype=float)
    ordered_labels = list('123456')
    selected_labels = list('123')
    # define the form objects
    form_objects = [
            Form.Matrix('matrix', 'distance matrix',
                D, MatrixUtil.assert_predistance),
            Form.MultiLine('labels', 'ordered labels',
                '\n'.join(ordered_labels)),
            Form.MultiLine('selection', 'selected labels',
                '\n'.join(selected_labels))]
    return form_objects

def set_to_string(my_set):
    """
    @param my_set: a sequence of arbitrary elements
    """
    return '{' + ', '.join(sorted(str(el) for el in my_set)) + '}'

def split_to_string(my_split):
    """
    @param my_split: a sequence of two taxon sequences
    """
    return ', '.join(set_to_string(taxa) for taxa in my_split)

def get_R_alpha_beta_gamma(D, index_set):
    """
    Four values are returned.
    R maps an index to the expected distance to a taxon in the other group;
    alpha gives the average distance from A to the defining branch;
    beta gives the average distance from B to the defining branch;
    gamma gives the length of the branch.
    @param D: a distance matrix
    @param index_set: the set of indices in group B
    @return: R, alpha, beta, gamma
    """
    n = len(D)
    A = set(range(n)) - set(index_set)
    B = set(index_set)
    nA = len(A)
    nB = len(B)
    if nA < 2 or nB < 2:
        raise ValueError('expected each side of the split to have at least two elements')
    R = {}
    R.update((i, sum(D[i,b] for b in B)/float(nB)) for i in A)
    R.update((j, sum(D[a,j] for a in A)/float(nA)) for j in B)
    gamma_plus_beta = 0.5 * min(R[i]+R[j]-D[i,j] for i, j in itertools.combinations(A, 2))
    alpha_plus_gamma = 0.5 * min(R[i]+R[j]-D[i,j] for i, j in itertools.combinations(B, 2))
    alpha_plus_gamma_plus_beta = sum(D[i,j] for i, j in itertools.product(A, B)) / float(nA * nB)
    gamma = alpha_plus_gamma + gamma_plus_beta - alpha_plus_gamma_plus_beta
    beta = gamma_plus_beta - gamma
    alpha = alpha_plus_gamma - gamma
    return R, alpha, beta, gamma

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # read the matrix
    D = np.array(fs.matrix)
    n = len(D)
    # read the ordered labels
    ordered_labels = Util.get_stripped_lines(StringIO(fs.labels))
    selected_labels = Util.get_stripped_lines(StringIO(fs.selection))
    # validate the input
    if n != len(ordered_labels):
        raise HandlingError('the number of taxon labels should match the number of rows in the distance matrix')
    # get the two sets of indices
    index_set_A = set(i for i, label in enumerate(ordered_labels) if label in selected_labels)
    index_set_B = set(range(n)) - index_set_A
    # get internal values related to the split
    R, alpha, beta, gamma = get_R_alpha_beta_gamma(D, index_set_B)
    # get the two new distance matrices
    D_A = BuildTreeTopology.update_generalized_nj(D, index_set_B)
    D_B = BuildTreeTopology.update_generalized_nj(D, index_set_A)
    # get the names associated with the indices of the new distance matrices
    all_names = [set([name]) for name in ordered_labels]
    D_A_names = [set_to_string(x) for x in SchurAlgebra.vmerge(all_names, index_set_B)]
    D_B_names = [set_to_string(x) for x in SchurAlgebra.vmerge(all_names, index_set_A)]
    # show the results
    out = StringIO()
    print >> out, 'alpha:', alpha
    print >> out, 'beta:', beta
    print >> out, 'gamma:', gamma
    print >> out
    print >> out, 'new distance matrix corresponding to the selected names:'
    print >> out, MatrixUtil.m_to_string(D_A)
    print >> out
    print >> out, 'ordered labels corresponding to this matrix:'
    for name in D_A_names:
        print >> out, name
    print >> out
    print >> out, 'new distance matrix corresponding to the non-selected names:'
    print >> out, MatrixUtil.m_to_string(D_B)
    print >> out
    print >> out, 'ordered labels corresponding to this matrix:'
    for name in D_B_names:
        print >> out, name
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
