"""Search for a counterexample to a conjecture about Euclidean distance matrices.

The conjecture is that when you convert an EDM to a Laplacian-like matrix,
and then merge a pair of rows and columns by summation,
then the eigenvalues of this matrix will remain non-negative.
One of the eigenvalues is guaranteed to remain zero,
and its corresponding eigenvector is guaranteed to be the unit vector.
"""

from StringIO import StringIO
import time
import random

import numpy
from numpy import linalg

from SnippetUtil import HandlingError
import MatrixUtil
import SchurAlgebra
import Euclid
import Form

def get_form():
    """
    @return: the body of a form
    """
    return [Form.Integer('npoints', 'use this many points per matrix', 8, low=3, high=20)]

def sample_points(npoints):
    """
    Each generated point is in N-1 dimensional space if there are N points.
    @param npoints: use this many points
    @return: a list of points in Euclidean space
    """
    mu = 0
    sigma = 1
    points = []
    for i in range(npoints):
        points.append([random.gauss(mu, sigma) for j in range(npoints-1)])
    return points

def points_to_edm(points):
    """
    @param points: a list of points in Euclidean space
    @return: the matrix of pairwise squared Euclidean distances
    """
    n = len(points)
    D = numpy.zeros((n,n))
    for i, a in enumerate(points):
        for j, b in enumerate(points):
            if i != j:
                sqdist = sum((y-x)**2 for x, y in zip(a, b))
                D[i][j] = sqdist
    return D


class Counterexample:
    def __init__(self, points, D, L_eigenvalues, D_small):
        self.points = points
        self.D = D
        self.L_eigenvalues = L_eigenvalues
        self.D_small = D_small


def process(npoints, nseconds):
    """
    @param npoints: attempt to form each counterexample from this many points
    @param nseconds: allow this many seconds to run
    @return: a multi-line string that summarizes the results
    """
    start_time = time.time()
    best_result = None
    nchecked = 0
    while time.time() - start_time < nseconds:
        # look for a counterexample
        points = sample_points(npoints)
        D = points_to_edm(points)
        L = Euclid.edm_to_laplacian(D)
        L_small = SchurAlgebra.mmerge(L, set([0, 1]))
        w = linalg.eigvalsh(L_small)
        D_small = Euclid.laplacian_to_edm(L_small)
        result = Counterexample(points, D, w, D_small)
        # see if the counterexample is interesting
        if best_result is None:
            best_result = result
        elif min(result.L_eigenvalues) < min(best_result.L_eigenvalues):
            best_result = result
        nchecked += 1
    out = StringIO()
    print >> out, 'checked', nchecked, 'matrices each formed from', npoints, 'points'
    print >> out
    print >> out, 'eigenvalues of the induced matrix with lowest eigenvalue:'
    for value in reversed(sorted(best_result.L_eigenvalues)):
        print >> out, value
    print >> out
    print >> out, 'corresponding induced distance matrix:'
    print >> out, MatrixUtil.m_to_string(best_result.D_small)
    print >> out
    print >> out, 'the original distance matrix corresponding to this matrix:'
    print >> out, MatrixUtil.m_to_string(best_result.D)
    print >> out
    print >> out, 'the points that formed the original distance matrix:'
    for point in best_result.points:
        print >> out, '\t'.join(str(x) for x in point)
    return out.getvalue().strip()

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    nseconds = 2
    npoints = fs.npoints
    response_text = process(npoints, nseconds)
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, response_text

def main():
    nseconds = 5
    npoints = 8
    print process(npoints, nseconds)

if __name__ == '__main__':
    main()
