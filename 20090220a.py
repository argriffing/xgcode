"""Given a distance matrix, get the first two principal coords of each point.

Given a distance matrix,
get the first two principal coordinates of each point.
If the tree-like option is used,
then the squared pairwise distances between the output points approximates
the input distances.
If the plane-like option is used,
then the pairwise distances between the output points approximates
the input distances.
"""

from StringIO import StringIO
import math

import numpy as np

from SnippetUtil import HandlingError
import MatrixUtil
import Form
import FormOut

def get_form():
    """
    @return: the body of a form
    """
    # define the default distance matrix
    D = np.array([
            [0, 405, 278, 502, 414],
            [405, 0, 542, 521, 774],
            [278, 542, 0, 246, 137],
            [502, 521, 246, 0, 304],
            [414, 774, 137, 304, 0]])
    # define the form objects
    form_objects = [
            Form.Matrix('matrix', 'distance matrix',
                D, MatrixUtil.assert_predistance),
            Form.RadioGroup('options', 'topology', [
                Form.RadioItem('treelike', 'distances are tree-like', True),
                Form.RadioItem('planelike', 'distances are plane-like')])]
    return form_objects

def get_form_out():
    return FormOut.Report()

#FIXME use eigutil

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # read the distance matrix
    D = fs.matrix
    # if the distances are plane-like then square each element of the distance matrix
    if fs.planelike:
        D = D**2
    # get the A matrix
    A = -0.5 * D
    # get the doubly centered A matrix
    HAH = MatrixUtil.double_centered(A)
    # do the eigendecomposition
    eigenvalues, eigenvector_transposes = np.linalg.eigh(HAH)
    eigenvectors = eigenvector_transposes.T
    eigensystem = [(abs(w), w, v.tolist()) for w, v in zip(eigenvalues, eigenvectors)]
    sorted_eigensystem = list(reversed(sorted(eigensystem)))
    sorted_abs_eigenvalues, sorted_eigenvalues, sorted_eigenvectors = zip(*sorted_eigensystem)
    # get the points that approximate the distance matrix
    top_eigenvalues = sorted_eigenvalues[:2]
    top_eigenvectors = sorted_eigenvectors[:2]
    for eigenvalue in top_eigenvalues:
        if eigenvalue <= 0:
            msg = 'one of the top two eigenvalues was non-positive'
            raise HandlingError(msg)
    axes = []
    for eigenvalue, eigenvector in zip(top_eigenvalues, top_eigenvectors):
        axis = [v * math.sqrt(eigenvalue) for v in eigenvector]
        axes.append(axis)
    points = zip(*axes)
    # begin the response
    out = StringIO()
    for point in points:
        print >> out, '\t'.join(str(v) for v in point)
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
