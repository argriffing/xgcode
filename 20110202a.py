"""
Check an asymptotic formula for eigenvectors of a path Laplacian.

Assume that the path has unit length.
"""

import math
from StringIO import StringIO

import numpy as np

import Form
import FormOut
import Euclid
import EigUtil
import iterutils

def get_form():
    """
    @return: a list of form objects
    """
    # define the form objects
    form_objects = [
            Form.Integer('nvertices', 'use this many vertices',
                20, low=2, high=200),
            Form.Integer('naxes', 'use this many eigenvectors',
                4, low=1)]
    return form_objects

def get_form_out():
    return FormOut.Report()

def create_laplacian_matrix(nvertices):
    """
    @param affinity: affinity between adjacent vertices
    @param nvertices: the number of vertices in the graph
    @return: a numpy matrix
    """
    affinity = float(nvertices - 1)
    A = np.zeros((nvertices, nvertices), dtype=float)
    for i, j in iterutils.pairwise(range(nvertices)):
        A[i,j] = affinity
        A[j,i] = affinity
    L = Euclid.adjacency_to_laplacian(A)
    return L

def sinusoidal_approximation_b(N, n, k):
    """
    This is a simplified formula for the sinusoidal approximation.
    @param N: the total number of vertices
    @param n: 1 for the most interesting eigenpair
    @param k: 0 for the first path vertex
    """
    return math.cos(n * k * math.pi / float(N - 1))

def sinusoidal_approximation(N, n, k):
    """
    @param N: the total number of vertices
    @param n: 1 for the most interesting eigenpair
    @param k: 0 for the first path vertex
    """
    nsegments = N-1
    segment_length = 1.0 / nsegments
    t = k * segment_length
    x = -math.pi / 2 + t * math.pi
    return math.cos(n * (x + 0.5 * math.pi))

def get_response_content(fs):
    # check input compatibility
    if fs.nvertices < fs.naxes+1:
        raise ValueError(
                'attempting to plot too many eigenvectors '
                'for the given number of vertices')
    # construct the path Laplacian matrix
    N = fs.nvertices
    L = create_laplacian_matrix(N)
    # compute the eigendecomposition
    ws, vs = EigUtil.eigh(L)
    # reorder the eigenvalues and eigenvectors
    ws = ws[:-1][::-1]
    vs = vs[:-1][::-1]
    # write the report
    np.set_printoptions(linewidth=200, threshold=10000)
    out = StringIO()
    for i in range(fs.naxes):
        w = ws[i]
        v = vs[i]
        n = i+1
        #scaled_eigenvector = v / math.sqrt(w)
        scaled_eigenvector = v * math.sqrt(N * 0.5)
        print >> out, scaled_eigenvector
        prediction = np.array([
            sinusoidal_approximation_b(N, n, k) for k in range(N)])
        print >> out, prediction
        print >> out, scaled_eigenvector / prediction
        print >> out
    return out.getvalue()
