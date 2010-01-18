"""Look for a counterexample to the eigenvector sign sufficiency conjecture.

Are the signs of eigenvector loadings of -(1/2)HDH sufficient
to get the topology of a bifurcating tree?
"""

import StringIO
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

g_epsilon = 1e-10

class CounterexampleError(Exception): pass

def get_form():
    """
    @return: a list of form objects
    """
    # define the list of form objects
    form_objects = [
            Form.Integer('ntaxa', 'number of taxa', 5, low=3, high=20)]
    return form_objects

def sample_branch_lengths(tree):
    """
    Modify the tree by setting branch lengths.
    @param tree: a tree
    """
    for node in tree.preorder():
        if not node.is_root():
            branch_length = float(random.randrange(1, 1000))
            node.set_branch_length(branch_length)

def reset_branch_lengths(tree):
    """
    Set each branch length to the unit length.
    @param tree: a tree
    """
    for node in tree.preorder():
        if not node.is_root():
            branch_length = 1
            node.set_branch_length(branch_length)

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

def point_to_orthant(p):
    """
    @param p: a point in general position in Euclidean space
    @return: a tuple of elements in {1, -1}
    """
    return tuple(1 if x > 0 else -1 for x in p)

def process(nseconds, ntaxa):
    """
    @param nseconds: allow this many seconds to run or None to run forever
    @return: a multi-line string that summarizes the results
    """
    start_time = time.time()
    nsamples_rejected = 0
    nsamples_accepted = 0
    pattern_to_topo_surrogate = {}
    pattern_to_tree_string = {}
    counterexample_message = 'no counterexample was found'
    try:
        while True:
            elapsed_time = time.time() - start_time
            if nseconds and elapsed_time > nseconds:
                break
            # sample an xtree topology
            xtree = TreeSampler.sample_agglomerated_tree(ntaxa)
            # convert the xtree to a FelTree, although I guess this might not be necessary
            tree_string = xtree.get_newick_string()
            tree = NewickIO.parse(tree_string, FelTree.NewickTree)
            # get ordered ids and the number of leaves and some auxiliary variables
            ordered_ids = get_ordered_ids(tree)
            nleaves = len(list(tree.gen_tips()))
            id_to_index = dict((myid, i) for i, myid in enumerate(ordered_ids))
            # force every branch length to be the unit length
            reset_branch_lengths(tree)
            # get the unweighted distance matrix among tips in convenient hashable form
            D_unit = np.array(tree.get_partial_distance_matrix(ordered_ids))
            topo_surrogate = tuple(tuple(row.tolist()) for row in D_unit)
            # sample random branch lengths
            sample_branch_lengths(tree)
            # get the weighted tree string
            weighted_tree_string = NewickIO.get_newick_string(tree)
            # get the distance matrix relating the leaves
            D = np.array(tree.get_partial_distance_matrix(ordered_ids))
            # get the projections onto the MDS axes of the leaves
            X = Euclid.edm_to_points(D)
            # if any coordinate is near zero then reject the sample
            if np.min(np.abs(X)) < g_epsilon:
                nsamples_rejected += 1
                continue
            # do an orthogonal transformation that puts the first point in the positive orthant
            canonizing_vector = np.array(point_to_orthant(X[0]))
            X *= canonizing_vector
            # get the canonical sign pattern
            sign_pattern = tuple(point_to_orthant(row) for row in X)
            # compare the topo surrogate of this sign pattern to the one in memory
            expected_topo_surrogate = pattern_to_topo_surrogate.get(sign_pattern, None)
            if expected_topo_surrogate:
                if topo_surrogate != expected_topo_surrogate:
                    remembered_tree_string = pattern_to_tree_string[sign_pattern]
                    msg = 'these trees have the same sign pattern but different topologies: {%s, %s}' % (weighted_tree_string, remembered_tree_string)
                    raise CounterexampleError(msg)
            else:
                pattern_to_topo_surrogate[sign_pattern] = topo_surrogate
                pattern_to_tree_string[sign_pattern] = weighted_tree_string
            # increment the count of accepted samples
            nsamples_accepted += 1
    except KeyboardInterrupt, e:
        pass
    except CounterexampleError, e:
        counterexample_message = str(e)
    # make the response
    out = StringIO.StringIO()
    print >> out, elapsed_time, 'seconds elapsed'
    print >> out, '%d-leaf trees were sampled' % ntaxa
    print >> out, 'epsilon for nonzero coordinates:', g_epsilon
    print >> out, 'search status:', counterexample_message
    print >> out, nsamples_rejected, 'samples were rejected because of a near-zero coordinate'
    print >> out, nsamples_accepted, 'samples were accepted'
    print >> out, len(pattern_to_topo_surrogate), 'unique canonical sign patterns were observed'
    return out.getvalue().strip()

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # on the web we have a short attention span
    nseconds = 2
    # get arguements from the user
    ntaxa = fs.ntaxa
    # get the response
    result_string = process(nseconds, ntaxa)
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
    print process(args.nseconds, args.ntaxa)

if __name__ == '__main__': 
    parser = argparse.ArgumentParser(description=SnippetUtil.docstring_to_title(__doc__))
    parser.add_argument('--nseconds', type=int, default=0, help='seconds to run or 0 to run until ctrl-c') 
    parser.add_argument('--ntaxa', type=int, default=5, help='number of taxa in each sampled tree topology') 
    args = parser.parse_args() 
    main(args) 
