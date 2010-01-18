"""Look for a counterexample to the generalized Fiedler eigenvector conjecture.

In a 2005 paper Fiedler looks at the hyperplane orthogonal to the nth principal
axis of the Steiner circumscribed hyperellipse bounding the simplex of a tree embedded in Euclidean space.
He notes that this hyperplane intersects the embedded tree at exactly n points.
In this snippet I look at whether this intersection property is true of the
hyperplanes of the Steiner circumscribed hyperellipse that bounds the simplex of
only the leaves of the tree embedded in Euclidean space.
"""

import StringIO
import math
import random
import time

import numpy as np

from SnippetUtil import HandlingError
import SnippetUtil
import Form
import MatrixUtil
import NewickIO
import FelTree
import Euclid
import TreeSampler
import TreeComparison


g_loading_epsilon = 1e-10


def get_form():
    """
    @return: a list of form objects
    """
    # define the list of form objects
    form_objects = []
    return form_objects

def do_internal_projection(D_full):
    """
    The resulting points are in the subspace whose basis vectors are the principal axes of the whole ellipsoid.
    @param D_full: the distance matrix as a numpy array relating all vertices including internal vertices
    @return: a numpy array where each row is a vertex of the tree
    """
    # Get the points such that the n rows in are points in n-1 dimensional space.
    # The first coordinate is the principal axis.
    points = Euclid.edm_to_points(D_full)
    return points

def do_external_projection(D_full, nleaves):
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

def process(nseconds):
    """
    @param nseconds: allow this many seconds to run or None to run forever
    @return: a multi-line string that summarizes the results
    """
    start_time = time.time()
    nsampled_trees = 0
    # track the number of observations of each number of cuts on each axis for each hyperellipse
    internal_important_axis_to_ncuts_dict = {}
    internal_unimportant_axis_to_ncuts_dict = {}
    external_axis_to_ncuts_dict = {}
    # track the number of bad axes of each principality for each hyperellipse
    internal_important_bad_axis_dict = {}
    internal_unimportant_bad_axis_dict = {}
    external_bad_axis_dict = {}
    try:
        while True:
            elapsed_time = time.time() - start_time
            if nseconds and elapsed_time > nseconds:
                break
            # pick a random number of taxa to use as leaves in the tree
            ntaxa = random.randrange(3, 12)
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
            # compute the set of pairs of indices corresponding to branches
            neighbor_index_pairs = set()
            for parent in tree.preorder():
                for child in parent.gen_children():
                    parent_index = id_to_index[id(parent)]
                    child_index = id_to_index[id(child)]
                    index_pair = frozenset((parent_index, child_index))
                    neighbor_index_pairs.add(index_pair)
            # get the distance matrix relating all of the points
            D_full = np.array(tree.get_full_distance_matrix(ordered_ids))
            # analyze the intersections of the axes of the ellipsoid that includes internal points
            internal_projection = do_internal_projection(D_full)
            npoints, naxes = internal_projection.shape
            # analyze low axes
            for axis in range(0, nleaves-1):
                if any(abs(internal_projection[i, axis]) < g_loading_epsilon for i in range(npoints)):
                    internal_important_bad_axis_dict[axis] = internal_important_bad_axis_dict.get(axis, 0) + 1
                else:
                    ncuts = 0
                    for indexa, indexb in neighbor_index_pairs:
                        if internal_projection[indexa, axis] * internal_projection[indexb, axis] < 0:
                            ncuts += 1
                    ncuts_dict = internal_important_axis_to_ncuts_dict.get(axis, {})
                    ncuts_dict[ncuts] = ncuts_dict.get(ncuts, 0) + 1
                    internal_important_axis_to_ncuts_dict[axis] = ncuts_dict
            # analyze high axes
            for axis in range(nleaves-1, naxes):
                if any(abs(internal_projection[i, axis]) < g_loading_epsilon for i in range(npoints)):
                    internal_unimportant_bad_axis_dict[axis] = internal_unimportant_bad_axis_dict.get(axis, 0) + 1
                else:
                    ncuts = 0
                    for indexa, indexb in neighbor_index_pairs:
                        if internal_projection[indexa, axis] * internal_projection[indexb, axis] < 0:
                            ncuts += 1
                    ncuts_dict = internal_unimportant_axis_to_ncuts_dict.get(axis, {})
                    ncuts_dict[ncuts] = ncuts_dict.get(ncuts, 0) + 1
                    internal_unimportant_axis_to_ncuts_dict[axis] = ncuts_dict
            # analyze the intersections of the axes of the ellipsoid that includes only leaf points
            external_projection = do_external_projection(D_full, nleaves)
            npoints, naxes = external_projection.shape
            for axis in range(naxes):
                if any(abs(external_projection[i, axis]) < g_loading_epsilon for i in range(npoints)):
                    external_bad_axis_dict[axis] = external_bad_axis_dict.get(axis, 0) + 1
                else:
                    ncuts = 0
                    for indexa, indexb in neighbor_index_pairs:
                        if external_projection[indexa, axis] * external_projection[indexb, axis] < 0:
                            ncuts += 1
                    ncuts_dict = external_axis_to_ncuts_dict.get(axis, {})
                    ncuts_dict[ncuts] = ncuts_dict.get(ncuts, 0) + 1
                    external_axis_to_ncuts_dict[axis] = ncuts_dict
            # increment the count of sampled trees
            nsampled_trees += 1
    except KeyboardInterrupt, e:
        pass
    # make the response
    out = StringIO.StringIO()
    print >> out, elapsed_time, 'seconds elapsed'
    print >> out, nsampled_trees, 'trees were sampled'
    print >> out, 'epsilon for ambiguous cuts:', g_loading_epsilon
    # show counts of ambiguous cuts
    print >> out, 'number of ambiguous important cuts by axis principality using the full ellipsoid:'
    if internal_important_bad_axis_dict:
        print >> out, ', '.join(str(axis+1) + ':' + str(count) for axis, count in sorted(internal_important_bad_axis_dict.items()))
    else:
        print >> out, '(none)'
    print >> out, 'number of ambiguous unimportant cuts by axis principality using the full ellipsoid:'
    if internal_unimportant_bad_axis_dict:
        print >> out, ', '.join(str(axis+1) + ':' + str(count) for axis, count in sorted(internal_unimportant_bad_axis_dict.items()))
    else:
        print >> out, '(none)'
    print >> out, 'number of ambiguous cuts by axis principality using the leaf ellipsoid:'
    if external_bad_axis_dict:
        print >> out, ', '.join(str(axis+1) + ':' + str(count) for axis, count in sorted(external_bad_axis_dict.items()))
    else:
        print >> out, '(none)'
    # show results for unambiguous cuts
    print >> out, 'number of observations of each number of cuts by each important axis of the full ellipsoid:'
    for axis, ncuts_dict in sorted(internal_important_axis_to_ncuts_dict.items()):
        print >> out, axis+1, ':', '{', ', '.join(str(ncuts) + ':' + str(count) for ncuts, count in sorted(ncuts_dict.items())), '}'
    print >> out, 'number of observations of each number of cuts by each unimportant axis of the full ellipsoid:'
    for axis, ncuts_dict in sorted(internal_unimportant_axis_to_ncuts_dict.items()):
        print >> out, axis+1, ':', '{', ', '.join(str(ncuts) + ':' + str(count) for ncuts, count in sorted(ncuts_dict.items())), '}'
    print >> out, 'number of observations of each number of cuts by each axis of the leaf ellipsoid:'
    for axis, ncuts_dict in sorted(external_axis_to_ncuts_dict.items()):
        print >> out, axis+1, ':', '{', ', '.join(str(ncuts) + ':' + str(count) for ncuts, count in sorted(ncuts_dict.items())), '}'
    return out.getvalue().strip()

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # on the web we have a short attention span
    nseconds = 2
    # get the response
    result_string = process(nseconds)
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

def main():
    nseconds = None
    print process(nseconds)

if __name__ == '__main__':
    main()
