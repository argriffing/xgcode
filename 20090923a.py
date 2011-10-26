"""Seek a counterexample to a principal orthant edge connectivity conjecture.

Look for a counterexample to the principal orthant
edge connectivity conjecture.
"""

from StringIO import StringIO
import random
import time

import numpy as np

from SnippetUtil import HandlingError
import SnippetUtil
import Form
import FormOut
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

def get_ncuts(pta, ptb):
    """
    @param pta: one endpoint in Euclidean space
    @param ptb: the other endpoint in Euclidean space
    """
    return sum(1 for a, b in zip(pta, ptb) if a*b < 0)

def get_internal_points(pta, ptb):
    """
    Get an ordered list of representative points in other orthants between two points.
    @param pta: the initial point in Euclidean space
    @param ptb: the final point in Euclidean space
    @return: an ordered list of points
    """
    # if there are fewer than two sign changes then nothing interesting is happening
    ncuts = get_ncuts(pta, ptb)
    if ncuts < 2:
        return []
    # if the midpoint is interesting then include it in the resulting list of points
    ptm = 0.5 * (pta + ptb)
    next_ncuts = get_ncuts(pta, ptm) 
    prefix = get_internal_points(pta, ptm)
    suffix = get_internal_points(ptm, ptb)
    if 0 < next_ncuts < ncuts:
        return prefix + [ptm] + suffix
    else:
        return prefix + suffix

def point_to_orthant(p):
    """
    @param p: a point in general position in Euclidean space
    @return: a tuple of elements in {1, -1}
    """
    return tuple(1 if x > 0 else -1 for x in p)

def get_blocked_orthants(pta, ptb):
    """
    Each blocked orthant contains an edge fragment.
    Each orthant is a tuple of elements in {1, -1}.
    @param pta: the initial point in Euclidean space
    @param ptb: the final point in Euclidean space
    @return: a set of orthants between pta and ptb exclusive
    """
    points = get_internal_points(pta, ptb)
    return set(point_to_orthant(p) for p in points)

def process(nseconds=None):
    """
    @param nseconds: allow this many seconds to run or None to run forever
    @return: a multi-line string that summarizes the results
    """
    start_time = time.time()
    nsampled_trees = 0
    counterexample_message = 'no counterexample was found'
    nvertex_connectivity_failures = 0
    nfragment_fragment_collisions = 0
    nfragment_vertex_collisions = 0
    ncontrol_vertex_connectivity_failures = 0
    ncontrol_fragment_fragment_collisions = 0
    ncontrol_fragment_vertex_collisions = 0
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
                npoints, naxes_total = projection.shape
                # look for a counterexample for each possible number of principal hyperplanes
                for naxes in range(1, naxes_total+1):
                    # some orthants are occupied by a fragment of an edge
                    forbidden_orthants = set()
                    for indexa, indexb in neighbor_index_pairs:
                        # get the endpoints of the edge in the Euclidean subspace
                        pta = projection[indexa][:naxes]
                        ptb = projection[indexb][:naxes]
                        # look at the orthants blocked by the fragments of this edge
                        orthants = get_blocked_orthants(pta, ptb)
                        if orthants & forbidden_orthants:
                            if D is D_full:
                                nfragment_fragment_collisions += 1
                                msg = 'two edge fragments occupy the same orthant in %d dimensions in the tree %s' % (naxes, tree_string)
                                raise CounterexampleError(msg)
                            else:
                                ncontrol_fragment_fragment_collisions += 1
                        forbidden_orthants.update(orthants)
                    # no vertex should share an orthant with an edge fragment
                    for i in range(npoints):
                        p = projection[i][:naxes]
                        orthant = point_to_orthant(p)
                        if orthant in forbidden_orthants:
                            if D is D_full:
                                nfragment_vertex_collisions += 1
                                msg = 'a vertex occupies the same orthant as an edge fragment in %d dimensions in the tree %s' % (naxes, tree_string)
                                raise CounterexampleError(msg)
                            else:
                                ncontrol_fragment_vertex_collisions += 1
                    # now partition the vertices by orthant and check their connectivity
                    orthant_to_id_set = {}
                    for i in range(npoints):
                        p = projection[i][:naxes]
                        orthant = point_to_orthant(p)
                        id_set = orthant_to_id_set.get(orthant, set())
                        id_set.add(ordered_ids[i])
                        orthant_to_id_set[orthant] = id_set
                    for id_set in orthant_to_id_set.values():
                        if not is_connected(tree, id_set):
                            if D is D_full:
                                nvertex_connectivity_failures += 1
                                msg = 'found disconnected vertices in an orthant in %d dimensions in the tree %s' % (naxes, tree_string)
                                raise CounterexampleError(msg)
                            else:
                                ncontrol_vertex_connectivity_failures += 1
            # increment the count of sampled trees
            nsampled_trees += 1
    except KeyboardInterrupt, e:
        pass
    except CounterexampleError as e
        counterexample_message = str(e)
    # make the response
    out = StringIO()
    print >> out, elapsed_time, 'seconds elapsed'
    print >> out, nsampled_trees, 'trees were sampled'
    print >> out, 'counterexample search status:', counterexample_message
    print >> out, 'vertex connectivity failures:', nvertex_connectivity_failures
    print >> out, 'fragment-fragment orthant collisions:', nfragment_fragment_collisions
    print >> out, 'fragment-vertex orthant collisions:', nfragment_vertex_collisions
    print >> out, 'control vertex connectivity failures:', ncontrol_vertex_connectivity_failures
    print >> out, 'control fragment-fragment orthant collisions', ncontrol_fragment_fragment_collisions
    print >> out, 'control fragment-vertex orthant collisions', ncontrol_fragment_vertex_collisions
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
