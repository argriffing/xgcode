"""Given a Laplacian matrix, find corresponding Euclidean points.

Given a Laplacian matrix, find points that define a Euclidean distance matrix.
Each row in the output defines a point in Euclidean space.
Contrasts are available for both normalized and combinatorial Laplacians,
while Euclidean coordinates are with respect to
the combinatorial Laplacian matrix.
"""

from StringIO import StringIO
import math

import numpy as np

from SnippetUtil import HandlingError
import MatrixUtil
from Form import RadioItem
from Form import CheckItem
import Form
import FormOut

def get_form():
    """
    @return: the body of a form
    """
    L = np.array([
            [1.0, 0.0, 0.0, 0.0, -1.0, 0.0],
            [0.0, 1.0, 0.0, 0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0, -1.0],
            [0.0, 0.0, 0.0, 1.0, 0.0, -1.0],
            [-1.0, -1.0, 0.0, 0.0, 3.0, -1.0],
            [0.0, 0.0, -1.0, -1.0, -1.0, 3.0]])
    # define the form objects
    form_objects = [
            Form.Matrix('laplacian', 'combinatorial Laplacian matrix',
                L, MatrixUtil.assert_symmetric),
            Form.RadioGroup('normalization', 'use this matrix', [
                RadioItem('combinatorial', 'combinatorial Laplacian', True),
                RadioItem('normalized', 'normalized Laplacian')]),
            Form.CheckGroup('options', 'output options', [
                CheckItem('show_magnitudes',
                    'show the distance of each point from the origin'),
                CheckItem('show_contrasts',
                    'show the contrast loading matrix')])]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_eigendecomposition(M):
    """
    @param M: a numpy array
    @return: the eigenvalues and the eigenvectors
    """
    w, v = np.linalg.eigh(M)
    eigenvalues = w
    eigenvectors = v.T
    return eigenvalues, eigenvectors

def gen_euclidean_points_from_eigendecomposition(w, v):
    """
    Yield Euclidean points from the eigendecomposition of the Laplacian.
    @param w: eigenvalues
    @param v: eigenvectors
    """
    epsilon = 1e-10
    if sum(1 for eigenvalue in w if abs(eigenvalue) < epsilon) != 1:
        msg = 'expected a single eigenvalue to be near zero: ' + str(w)
        raise ValueError(msg)
    n = len(w)
    for i in range(n):
        p = []
        for j in range(n):
            eigenvalue = w[j]
            if abs(eigenvalue) > epsilon:
                coord = v[j][i] / math.sqrt(eigenvalue)
                p.append(coord)
        yield p

def gen_contrasts_from_eigendecomposition(w, v):
    """
    Yield contrast loading vectors from the eigendecomposition of the Laplacian.
    @param w: eigenvalues
    @param v: eigenvectors
    """
    epsilon = 1e-10
    if sum(1 for eigenvalue in w if abs(eigenvalue) < epsilon) != 1:
        raise ValueError('expected a single eigenvalue to be near zero')
    n = len(w)
    for i in range(n):
        p = []
        for j in range(n):
            eigenvalue = w[j]
            if abs(eigenvalue) > epsilon:
                coord = v[j][i] * math.sqrt(eigenvalue)
                p.append(coord)
        yield p

def gen_normalized_laplacian_contrasts(w, v, d):
    """
    Yield contrast loading vectors from the eigendecomposition of the normalized Laplacian.
    @param w: eigenvalues
    @param v: eigenvectors
    @param d: the row sums of the affinity matrix
    """
    epsilon = 1e-10
    if sum(1 for eigenvalue in w if abs(eigenvalue) < epsilon) != 1:
        raise ValueError('expected a single eigenvalue to be near zero')
    n = len(w)
    for i in range(n):
        degree = d[i]
        p = []
        for j in range(n):
            eigenvalue = w[j]
            if abs(eigenvalue) > epsilon:
                coord = math.sqrt(degree) * v[j][i] * math.sqrt(eigenvalue)
                p.append(coord)
        yield p


def get_response_content(fs):
    # read the laplacian matrix
    L = fs.laplacian
    n = len(L)
    # get the eigenvalues and the eigenvectors of the Laplacian
    eigenvalues, eigenvectors = get_eigendecomposition(L)
    # get the Euclidean points
    points = list(gen_euclidean_points_from_eigendecomposition(eigenvalues, eigenvectors))
    distances = [sum(v*v for v in point) for point in points]
    # get the contrasts from the combinatorial or from the normalized laplacian
    if fs.combinatorial:
        contrasts = list(gen_contrasts_from_eigendecomposition(eigenvalues, eigenvectors))
    else:
        # get the row sums of the affinity matrix
        row_sums = [L[i][i] for i in range(n)]
        # get the normalized laplacian
        L_script = np.zeros((n,n))
        for i in range(n):
            for j in range(n):
                L_script[i][j] = L[i][j] / math.sqrt(row_sums[i]*row_sums[j])
        # get the eigenvalues and the eigenvectors of the normalized Laplacian
        eigenvalues, eigenvectors = get_eigendecomposition(L_script)
        # get the contrasts
        contrasts = list(gen_normalized_laplacian_contrasts(eigenvalues, eigenvectors, row_sums))
    # define the paragraphs of the response
    paragraphs = []
    if True:
        paragraph = []
        if fs.normalized:
            paragraph.append('warning: using the combinatorial laplacian for the points')
        paragraph.extend([
            'each row is a point:',
            '\n'.join('\t'.join(str(value) for value in point) for point in points)])
        paragraphs.append(paragraph)
    if fs.show_magnitudes:
        paragraph = [
                'distance from each point to the origin:',
                '\n'.join(str(d) for d in distances)]
        paragraphs.append(paragraph)
    if fs.show_contrasts:
        paragraph = [
                'contrast loadings:',
                '\n'.join('\t'.join(str(value) for value in point) for point in contrasts)]
        paragraphs.append(paragraph)
    # return the response
    return '\n\n'.join('\n'.join(p) for p in paragraphs) + '\n'
