"""Look for a counterexample to a principal orthant connectivity conjecture.

Look for a counterexample to the principal orthant
vertex connectivity conjecture.
"""

from StringIO import StringIO
import random
import time

import numpy as np

from SnippetUtil import HandlingError
import SnippetUtil
import Form
import FormOut
import MatrixUtil
import NewickIO
import FelTree
import Euclid
import TreeSampler
import TreeComparison

def CounterexampleError(Exception): pass

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

def is_connected(tree, ids):
    """
    This uses an iterative breadth first enumeration.
    @param tree: a tree
    @param ids: a set of ids of nodes whose connectivity in the tree is being queried
    @return: True or False depending on whether the ids are connected in the tree
    """
    # get the mapping from id to node
    id_to_node = dict((id(node), node) for node in tree.preorder())
    # initialize the shell with an arbitrary id in the set
    next_shell = set([list(ids)[0]])
    observed = set(next_shell)
    while next_shell:
        shell = next_shell
        next_shell = set()
        for myid in shell:
            for neighbor in id_to_node[myid].gen_neighbors():
                nextid = id(neighbor)
                if nextid in ids and nextid not in observed:
                    observed.add(nextid)
                    next_shell.add(nextid)
    return (observed == ids)

def process(nseconds=None):
    """
    @param nseconds: allow this many seconds to run or None to run forever
    @return: a multi-line string that summarizes the results
    """
    start_time = time.time()
    nsampled_trees = 0
    counterexample_message = 'no counterexample was found'
    northants_passed = 0
    northants_failed = 0
    ncontrol_orthants_passed = 0
    ncontrol_orthants_failed = 0
    branch_cut_hist = {}
    control_branch_cut_hist = {}
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
            D_bad = np.sqrt(D_full)
            for D in (D_full, D_bad):
                # get the projections onto the MDS axes of the leaves
                projection = do_projection(D, nleaves)
                npoints, naxes = projection.shape
                # recursively split the points by hyperplanes of principal axes
                next_id_set_list = [set(ordered_ids)]
                for axis in range(naxes):
                    id_set_list = next_id_set_list
                    # create the list of sets of points in principal orthants
                    next_id_set_list = []
                    for id_set in id_set_list:
                        neg_id_set = set(myid for myid in id_set if projection[id_to_index[myid], axis] < 0)
                        nonneg_id_set = set(myid for myid in id_set if projection[id_to_index[myid], axis] >= 0)
                        for next_set in (neg_id_set, nonneg_id_set):
                            if len(next_set) > 1:
                                next_id_set_list.append(next_set)
                    # each set of points should be connected
                    for id_set in next_id_set_list:
                        bconnected = is_connected(tree, id_set)
                        if bconnected and (D is D_full):
                            northants_passed += 1
                        elif (not bconnected) and (D is D_full):
                            northants_failed += 1
                            msg = 'found a counterexample in principal orthant %d of the tree %s' % (axis+1, tree_string)
                            raise CounterexampleError(msg)
                        elif bconnected and (D is not D_full):
                            ncontrol_orthants_passed += 1
                        elif (not bconnected) and (D is not D_full):
                            ncontrol_orthants_failed += 1
                # define the applicable histogram
                hist = branch_cut_hist if D is D_full else control_branch_cut_hist
                # check the number of cuts per branch
                for i, j in neighbor_index_pairs:
                    ncuts = sum(1 for axis in range(naxes) if projection[i, axis] * projection[j, axis] < 0)
                    hist[ncuts] = hist.get(ncuts, 0) + 1
            # increment the count of sampled trees
            nsampled_trees += 1
    except KeyboardInterrupt, e:
        pass
    except CounterexampleError as e:
        counterexample_message = str(e)
    # make the response
    out = StringIO()
    print >> out, elapsed_time, 'seconds elapsed'
    print >> out, nsampled_trees, 'trees were sampled'
    print >> out, 'counterexample search status:', counterexample_message
    print >> out, 'number of orthants with a connected component:', northants_passed
    print >> out, 'number of orthants with disconnected components:', northants_failed
    print >> out, 'number of control orthants with a connected component:', ncontrol_orthants_passed
    print >> out, 'number of control orthants with disconnected components:', ncontrol_orthants_failed
    print >> out, 'histogram of number of cuts per branch (real):'
    print >> out, ', '.join(str(ncuts) + ':' + str(count) for ncuts, count in branch_cut_hist.items())
    print >> out, 'histogram of number of cuts per branch (control):'
    print >> out, ', '.join(str(ncuts) + ':' + str(count) for ncuts, count in control_branch_cut_hist.items())
    return out.getvalue().strip()

def get_response_content(fs):
    # on the web we have a short attention span
    return process(nseconds=2) + '\n'

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
    print process()

if __name__ == '__main__':
    main()
