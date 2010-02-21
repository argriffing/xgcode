"""Given a tree, derive a rotated contrast matrix from its distance matrix.

Given a newick tree, derive a rotated contrast matrix from its distance matrix.
Elements in the output matrix with absolute values smaller than epsilon
will be zeroed.
A command in R to orthogonally rotate matrix M is
"varimax(M, normalize=FALSE)".
A command in MATLAB to orthogonally rotate matrix M is
"rotatefactors(M, 'method', 'varimax', 'normalize', 'off')".
"""

from StringIO import StringIO
import math

import numpy as np

from SnippetUtil import HandlingError
import Util
import NewickIO
import FelTree
import MatrixUtil
import Euclid
import Form

def get_form():
    """
    @return: the body of a form
    """
    # define the default tree string
    tree_string = '((((a:0.05, b:0.05):0.15, c:0.2):0.8, x:1.0):0.5, (((m:0.05, n:0.05):0.15, p:0.2):0.8, y:1.0):0.5);'
    tree = NewickIO.parse(tree_string, FelTree.NewickTree)
    formatted_tree_string = NewickIO.get_narrow_newick_string(tree, 60)
    # define the ordered labels
    ordered_labels = list('abcxmnpy')
    # define the form objects
    form_objects = [
            Form.MultiLine('tree', 'newick tree with branch lengths',
                formatted_tree_string),
            Form.MultiLine('labels', 'ordered labels',
                '\n'.join(ordered_labels)),
            Form.Float('epsilon', 'epsilon', '1e-10'),
            Form.RadioGroup('matrix_format', 'output matrix format', [
                Form.RadioItem('plain_format', 'plain', True),
                Form.RadioItem('r_format', 'R'),
                Form.RadioItem('matlab_format', 'MATLAB')])]
    return form_objects

def get_eigendecomposition(M):
    """
    @param M: a numpy array
    @return: the eigenvalues and the eigenvectors
    """
    w, v = np.linalg.eigh(M)
    eigenvalues = w
    eigenvectors = v.T
    return eigenvalues, eigenvectors

def get_contrast_matrix(w, v, eps=1e-10):
    """
    @param w: eigenvalues
    @param v: eigenvectors
    @param eps: exactly one eigenvalue should be closer than this to zero
    @return: an N x N-1 contrast matrix as a numpy array
    """
    abs_eigenvalues = [abs(eigenvalue) for eigenvalue in w]
    nsmall = sum(1 for abs_eigenvalue in abs_eigenvalues if abs_eigenvalue < eps)
    if not nsmall:
        raise ValueError('the minimum absolute value of any eigenvector is ' + str(min(abs_eigenvalues)))
    elif nsmall > 1:
        raise ValueError('observed %d eigenvalues near zero' % nsmall)
    n = len(w)
    points = []
    for i in range(n):
        p = []
        for j in range(n):
            eigenvalue = w[j]
            if abs(eigenvalue) > eps:
                coord = v[j][i] * math.sqrt(eigenvalue)
                p.append(coord)
        points.append(p)
    return np.array(points)

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # get the tree
    tree = NewickIO.parse(fs.tree, FelTree.NewickTree)
    # read the ordered labels
    ordered_labels = Util.get_stripped_lines(StringIO(fs.labels))
    # validate the input
    observed_label_set = set(node.get_name() for node in tree.gen_tips())
    if set(ordered_labels) != observed_label_set:
        raise HandlingError('the ordered labels should match the labels of the leaves of the tree')
    # get the matrix of pairwise distances among the tips
    D = np.array(tree.get_distance_matrix(ordered_labels))
    L = Euclid.edm_to_laplacian(D)
    w, v = get_eigendecomposition(L)
    C = get_contrast_matrix(w, v)
    # set elements with small absolute value to zero
    C[abs(C) < fs.epsilon] = 0
    # start to prepare the reponse
    out = StringIO()
    if fs.plain_format:
        print >> out, MatrixUtil.m_to_string(C)
    elif fs.matlab_format:
        print >> out, MatrixUtil.m_to_matlab_string(C)
    elif fs.r_format:
        print >> out, MatrixUtil.m_to_R_string(C)
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
