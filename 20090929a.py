"""Look at each step of the projection onto the axes of the Steiner ellipsoid of the leaves.

Because this really is a projection it can obviously be done in one step,
but I haven't worked this out yet.
"""

import StringIO
import random
import time
import argparse

import numpy as np

from SnippetUtil import HandlingError
import SnippetUtil
import Form
import NewickIO
import FelTree
import Euclid
import TreeSampler

def CounterexampleError(Exception): pass

def get_form():
    """
    @return: a list of form objects
    """
    # define the list of form objects
    form_objects = [
            Form.Integer('ntaxa', 'number of taxa', 5, low=3, high=20)]
    return form_objects

def process(ntaxa):
    np.set_printoptions(linewidth=200)
    out = StringIO.StringIO()
    # sample an xtree topology
    xtree = TreeSampler.sample_agglomerated_tree(ntaxa)
    # sample an xtree with exponentially distributed branch lengths
    mu = 2.0
    for branch in xtree.get_branches():
        branch.length = random.expovariate(1/mu)
    # convert the xtree to a FelTree so we can use the internal vertices
    tree_string = xtree.get_newick_string()
    tree = NewickIO.parse(tree_string, FelTree.NewickTree)
    # get ordered ids and the number of leaves and some auxiliary variables
    ordered_ids = get_ordered_ids(tree)
    nleaves = len(list(tree.gen_tips()))
    id_to_index = dict((myid, i) for i, myid in enumerate(ordered_ids))
    # get the distance matrix relating all of the points
    D_full = np.array(tree.get_full_distance_matrix(ordered_ids))
    # Now do the projection so that
    # the resulting points are in the subspace whose basis vectors are the axes of the leaf ellipsoid.
    # First get the points such that the n rows in X are points in n-1 dimensional space.
    X = Euclid.edm_to_points(D_full)
    print >> out, 'points with centroid at origin:'
    print >> out, X
    print >> out
    # Translate all of the points so that the origin is at the centroid of the leaves.
    X -= np.mean(X[:nleaves], 0)
    print >> out, 'points with centroid of leaves at origin:'
    print >> out, X
    print >> out
    # Extract the subset of points that define the leaves.
    L = X[:nleaves]
    # Find the orthogonal transformation of the leaves onto their MDS axes.
    # According to the python svd documentation, singular values are sorted most important to least important.
    U, s, Vt = np.linalg.svd(L)
    # Transform all of the points (including the internal vertices) according to this orthogonal transformation.
    # The axes are now the axes of the Steiner circumscribed ellipsoid of the leaf vertices.
    # I am using M.T[:k].T to get the first k columns of M.
    Z = np.dot(X, Vt.T)
    print >> out, 'orthogonally transformed points (call this Z):'
    print >> out, Z
    print >> out
    Y = Z.T[:(nleaves-1)].T
    print >> out, 'projection of the points onto the axes of the leaf ellipsoid,'
    print >> out, '(these are the first columns of Z; call this projected matrix Y):'
    print >> out, Y
    print >> out
    # Show the inner products.
    inner_products_of_columns = np.dot(Y.T, Y)
    print >> out, "pairwise inner products of the columns of Y (that is, Y'Y)"
    print >> out, inner_products_of_columns
    print >> out
    # Show other inner products.
    inner_products_of_columns = np.dot(Y[:5].T, Y[:5])
    print >> out, "pairwise inner products of the first few columns of Y"
    print >> out, inner_products_of_columns
    print >> out
    # Extract the subset of points that define the points of articulation.
    # Note that the origin is the centroid of the leaves.
    R = X[nleaves:]
    Y_leaves = Y[:nleaves]
    W = np.dot(np.linalg.pinv(L), Y_leaves)
    print >> out, 'leaf projection using pseudoinverse (first few rows of Y):'
    print >> out, np.dot(L, W)
    print >> out
    print >> out, 'projection of points of articulation using pseudoinverse (remaining rows of Y):'
    print >> out, np.dot(R, W)
    print >> out
    # Get all of the points in high dimensional space.
    X = Euclid.edm_to_points(D_full)
    # Get the MDS onto the lower dimensional space.
    X = X.T[:(nleaves-1)].T
    assert np.allclose(sum(X, 0), 0)
    print >> out, 'all points projected onto the first principal axes of the full ellipsoid:'
    print >> out, X
    print >> out
    # Look at only the leaves in this space.
    L = X[:nleaves]
    L -= np.mean(L, 0)
    print >> out, 'leaves projected onto the first principal axes of the full ellipsoid and then centered:'
    print >> out, L
    print >> out
    # Re-project the leaves onto the axes of leaf ellipsoid.
    D_leaves = Euclid.dccov_to_edm(np.dot(L, L.T))
    Y = Euclid.edm_to_points(D_leaves)
    print >> out, 'leaves further projected onto principal axes of their own ellipsoid:'
    print >> out, Y
    print >> out
    # Try something else
    D_all = Euclid.dccov_to_edm(np.dot(X, X.T))
    Y = Euclid.edm_to_points(D_all).T[:(nleaves-1)].T
    print >> out, 'all points further projected onto their own principal axes of inertia:'
    print >> out, Y
    print >> out
    # Try the same thing some more
    D_again = Euclid.dccov_to_edm(np.dot(Y, Y.T))
    Z = Euclid.edm_to_points(D_again).T[:(nleaves-1)].T
    print >> out, 'all points further projected onto their own principal axes of inertia (second iteration):'
    print >> out, Z
    print >> out
    return out.getvalue().strip()

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # get the response
    result_string = process(fs.ntaxa)
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, result_string

def get_ordered_ids(tree):
    """
    Maybe I could use postorder here instead.
    @param tree: a tree
    @return: a list of ids beginning with the leaves
    """
    ordered_ids = []
    ordered_ids.extend(id(node) for node in tree.gen_tips())
    ordered_ids.extend(id(node) for node in tree.gen_internal_nodes())
    return ordered_ids

def main(args):
    print process(args.ntaxa)

if __name__ == '__main__': 
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ntaxa', type=int, default=5, help='number of taxa in each sampled tree topology') 
    args = parser.parse_args() 
    main(args) 
