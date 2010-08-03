"""Compare two Euclidean embeddings of tree vertices.

The first embedding is the limit of the MDS coordinates
for all vertices in the tree as
masses of vertices of articulation go to zero
while masses of the leaves remain uniform.
The second embedding is the projection of the vertices of articulation
onto the MDS space defined by the distance matrix
of only the (uniformly weighted) leaves.
Are these two embeddings the same?
"""


from StringIO import StringIO
import random
import time

import numpy as np
import argparse

from SnippetUtil import HandlingError
import SnippetUtil
import Form
import FormOut
import NewickIO
import FelTree
import Euclid
import TreeSampler


def get_form():
    """
    @return: a list of form objects
    """
    # define the list of form objects
    form_objects = []
    return form_objects

def get_form_out():
    return FormOut.Report()

def do_projection(D_full, nleaves):
    """
    Project points onto the space of the leaves.
    The resulting points are in the subspace
    whose basis vectors are the principal axes of the leaf ellipsoid.
    @param D_full: distances relating all, including internal, vertices.
    @param nleaves: the first few indices in D_full represent leaves
    @return: a numpy array where each row is a vertex of the tree
    """
    # Get the points
    # such that the n rows in X are points in n-1 dimensional space.
    X = Euclid.edm_to_points(D_full)
    # Translate all of the points
    # so that the origin is at the centroid of the leaves.
    X -= np.mean(X[:nleaves], 0)
    # Extract the subset of points that define the leaves.
    L = X[:nleaves]
    # Find the orthogonal transformation of the leaves onto their MDS axes.
    # According to the python svd documentation,
    # singular values are sorted most important to least important.
    U, s, Vt = np.linalg.svd(L)
    # Transform all of the points (including the internal vertices)
    # according to this orthogonal transformation.
    # The axes are now the principal axes
    # of the Steiner circumscribed ellipsoid of the leaf vertices.
    # I am using M.T[:k].T to get the first k columns of M.
    points = np.dot(X, Vt.T).T[:(nleaves-1)].T
    return points

def process():
    """
    @return: a multi-line string that summarizes the results
    """
    np.set_printoptions(linewidth=200)
    out = StringIO()
    # define a degenerate mass vector
    m_degenerate = np.array([0.25, 0.25, 0.25, 0.25, 0, 0])
    # define some distance matrices
    D_leaves = Euclid.g_D_b
    D_all = Euclid.g_D_c
    nvertices = 6
    nleaves = 4
    # get the projection and the weighted multidimensional scaling
    X = Euclid.edm_to_points(D_all)
    Y = Euclid.edm_to_weighted_points(D_all, m_degenerate)
    D_X = np.array([[np.dot(pb-pa, pb-pa) for pa in X] for pb in X])
    D_Y = np.array([[np.dot(pb-pa, pb-pa) for pa in Y] for pb in Y])
    # get the embedding using only the leaves
    print >> out, 'embedding of leaves from the leaf distance matrix:'
    print >> out, Euclid.edm_to_points(D_leaves)
    print >> out, 'projection of all vertices onto the MDS space of the leaves:'
    print >> out, do_projection(D_all, nleaves)
    print >> out, 'embedding of all vertices using uniform weights:'
    print >> out, X
    print >> out, 'corresponding distance matrix:'
    print >> out, D_X
    print >> out, 'embedding of all vertices using degenerate weights:'
    print >> out, Y
    print >> out, 'corresponding distance matrix:'
    print >> out, D_Y
    return out.getvalue().strip()

def get_response_content(fs):
    return process() + '\n'

if __name__ == '__main__': 
    print process()
