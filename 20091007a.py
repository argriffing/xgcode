"""Look at the effect of added leaflets on the MDS embedding of a tree.

By a leaflet I mean a new leaf that is a negligible distance from an existing vertex.
Consider the following three situations:
1) No leaflet is added, in which case the embedding of the tree is not affected.
2) One leaflet is added to each leaf in the original tree.
3) A large number of leaflets are added to each leaf in the original tree.
I hope that case (3) is like computing the MDS using only leaf distances,
and then projecting the internal vertices onto this space.
Also I hope that case (2) can show that these leaflets can act as masses
that can be generalized to any real vector of homogeneous coordinates.
I had hoped that this generalization would have been given by the formulas
in the 2007 Abdi MDS chapter, but I think that this is not the case,
so now I will try to find the correct way to use the mass vector.
"""


from StringIO import StringIO
import random
import time

import numpy as np
import argparse

from SnippetUtil import HandlingError
import SnippetUtil
import Form
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

def add_leaflets(D_full, nleaves):
    """
    Create a modified distance matrix by effectively duplicating each leaf.
    @param D_full: the distance matrix as a numpy array relating all vertices including internal vertices
    @param nleaves: the first few indices in D_full represent leaves
    @return: a higher order distance matrix with new leaflet indices
    """
    # define the total number of vertices in the orginal tree
    nvertices = len(D_full)
    # define the total number of vertices in the next tree
    nnext = nvertices + nleaves
    # begin creating the next distance matrix
    D_next = np.zeros((nnext, nnext))
    # define distances among the original vertices
    for i in range(nvertices):
        for j in range(nvertices):
            D_next[i,j] = D_full[i,j]
    # define distances between old vertices and leaflets
    for i in range(nleaves):
        for j in range(nvertices):
            D_next[nvertices+i,j] = D_full[i,j]
            D_next[j,nvertices+i] = D_full[i,j]
    # define distances among leaflets
    for i in range(nleaves):
        for j in range(nleaves):
            D_next[nvertices+i,nvertices+j] = D_full[i,j]
    return D_next

def do_projection(D_full, nleaves):
    """
    The resulting points are in the subspace whose basis vectors are the principal axes of the leaf ellipsoid.
    @param D_full: the distance matrix as a numpy array relating all vertices including internal vertices
    @param nleaves: the first few indices in D_full represent leaves
    @return: a numpy array where each row is a vertex of the tree
    """
    # Get the points such that the n rows in X are points in n-1 dimensional space.
    X = Euclid.edm_to_points(D_full)
    # Translate all of the points so that the origin is at the centroid of the leaves.
    X -= np.mean(X[:nleaves], 0)
    # Extract the subset of points that define the leaves.
    L = X[:nleaves]
    # Find the orthogonal transformation of the leaves onto their MDS axes.
    # According to the python svd documentation, singular values are sorted most important to least important.
    U, s, Vt = np.linalg.svd(L)
    # Transform all of the points (including the internal vertices) according to this orthogonal transformation.
    # The axes are now the principal axes of the Steiner circumscribed ellipsoid of the leaf vertices.
    # I am using M.T[:k].T to get the first k columns of M.
    points = np.dot(X, Vt.T).T[:(nleaves-1)].T
    return points

def get_weighted_embedding_b(D, m):
    """
    This method was suggested by Eric but it has some limitations.
    In particular, it fails when some elements of the mass vector are zero.
    @param D: a distance matrix
    @param m: a mass vector
    @return: an embedding
    """
    M = np.diag(np.sqrt(m))
    cross_product_matrix = Euclid.edm_to_weighted_cross_product(D, m)
    Q = np.dot(M, np.dot(cross_product_matrix, M.T))
    U, S, VT = np.linalg.svd(Q, full_matrices=False)
    Z = np.dot(np.linalg.pinv(M), np.dot(U, np.sqrt(np.diag(S))))
    return Z

def process():
    """
    @return: a multi-line string that summarizes the results
    """
    np.set_printoptions(linewidth=200)
    out = StringIO()
    # define some distance matrices
    D_leaves = Euclid.g_D_b
    D_all = Euclid.g_D_c
    nvertices = 6
    nleaves = 4
    # define mass vectors
    m_degenerate = np.array([0.25, 0.25, 0.25, 0.25, 0, 0])
    m_interesting = np.array([.2, .2, .2, .2, .1, .1])
    m_uniform = np.ones(nvertices) / float(nvertices)
    # augment a distance matrix by adding leaflets
    D_augmented = add_leaflets(D_all, nleaves)
    # create the projection of points
    X_projected = do_projection(D_all, nleaves)
    # show some of the distance matrices
    print >> out, 'pairwise distances among vertices in the original tree:'
    print >> out, D_all
    print >> out, 'pairwise distance matrix augmented with one leaflet per leaf:'
    print >> out, D_augmented
    # get the distance matrices corresponding to the cases in the docstring
    print >> out, 'case 1: embedding of all vertices:'
    print >> out, Euclid.edm_to_points(D_all)
    print >> out, 'case 2: embedding of leaves and leaflets from the leaflet-augmented distance matrix:'
    print >> out, Euclid.edm_to_points(D_augmented)
    print >> out, 'case 3: projection of all vertices onto the MDS space of the leaves:'
    print >> out, X_projected
    # another embedding
    print >> out, 'embedding of leaves from the leaf distance matrix:'
    print >> out, Euclid.edm_to_points(D_leaves)
    # show embeddings of a tree augmented with leaflets
    print >> out, 'first few coordinates of the original vertices of the embedded tree with lots of leaflets per leaf:'
    D_super_augmented = D_all.copy()
    for i in range(20):
        D_super_augmented = add_leaflets(D_super_augmented, nleaves)
    X_super = Euclid.edm_to_points(D_super_augmented)
    X_super_block_small = X_super[:6].T[:3].T
    print >> out, X_super_block_small
    print >> out, 'ratio of coordinates of projected points to coordinates of this block of the embedding of the augmented tree:'
    print >> out, X_projected / X_super_block_small
    # test
    Z = Euclid.edm_to_weighted_points(D_all, m_uniform)
    print >> out, 'generalized case 1:'
    print >> out, Z
    # test
    Z = Euclid.edm_to_weighted_points(D_all, m_interesting)
    print >> out, 'generalized case 2:'
    print >> out, Z
    # test
    Z = Euclid.edm_to_weighted_points(D_all, m_degenerate)
    print >> out, 'generalized case 3:'
    print >> out, Z
    # test
    Z = get_weighted_embedding_b(D_all, m_uniform)
    print >> out, 'eric formula case 1:'
    print >> out, Z
    # test
    Z = get_weighted_embedding_b(D_all, m_interesting)
    print >> out, 'eric formula case 2:'
    print >> out, Z
    # test
    Z = get_weighted_embedding_b(D_all, m_degenerate)
    print >> out, 'eric formula case 3:'
    print >> out, Z
    # test stuff
    print >> out, 'testing random stuff:'
    D = D_all
    m = m_degenerate
    nvertices = len(m)
    sqrtm = np.sqrt(m)
    M = np.diag(sqrtm)
    cross_product_matrix = Euclid.edm_to_weighted_cross_product(D, m)
    U_cross, S_cross, VT_cross = np.linalg.svd(cross_product_matrix, full_matrices=False)
    Q = np.dot(M, np.dot(cross_product_matrix, M.T))
    U, B, VT = np.linalg.svd(Q, full_matrices=False)
    S = np.sqrt(np.diag(B))
    US = np.dot(U, S)
    M_pinv = np.linalg.pinv(M)
    M_pinv_narrow = M_pinv.T[:-2].T
    US_short = US[:-2]
    print >> out, 'eigenvalues of the abdi cross product:', S_cross
    print >> out, 'eigenvalues of the eric cross product:', B
    print >> out, M_pinv
    print >> out, US
    print >> out, M_pinv_narrow
    print >> out, US_short
    Z = np.dot(M_pinv_narrow, US_short)
    print >> out, Z
    # return the response
    return out.getvalue().strip()

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # get the response
    result_string = process()
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, result_string

if __name__ == '__main__': 
    print process()
