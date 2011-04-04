"""
Do an exhaustive search of internal vertex sign assignments.
"""

import unittest
from collections import defaultdict
import itertools

import iterutils


def is_branch_compat(nsame, ndifferent, ntarget, nbranches):
    """
    @param nsame: number of placed edges without sign change
    @param ndifferent: number of placed edges with sign change
    @param ntarget: target number of edges with sign change
    @param nbranches: the total number of branches in the tree.
    """
    if nsame + ndifferent > nbranches:
        raise ValueError('branch sign change error')
    if ndifferent > ntarget:
        return False
    npotential = nbranches - nsame
    if npotential < ntarget:
        return False
    return True

def rec_internal(
        id_to_adj, id_to_val,
        nsame, ndifferent, ntarget, nbranches,
        internals, depth):
    """
    This is a recursive function.
    Each level corresponds to an internal vertex.
    Each time a +1/-1 is assigned to an internal vertex,
    check that the number of sign changes on edges is correct.
    @param id_to_adj: node id to list of ids of adjacent nodes
    @param id_to_val: node id to valuation
    @param nsame: number of placed edges without sign change
    @param ndifferent: number of placed edges with sign change
    @param ntarget: target number of edges with sign change
    @param nbranches: the total number of branches in the tree.
    @param internals: list of ids of internal nodes
    @param depth: recursion depth starts at zero
    """
    idcur = internals[depth]
    for value in (-1, 1):
        # Check the number of edges where the signs
        # change and where the signs stay the same
        # under the proposed valuation for the current internal vertex.
        nsame_next = nsame
        ndifferent_next = ndifferent
        for adj in id_to_adj[idcur]:
            adj_val = id_to_val[adj]
            if adj_val is not None:
                prod = adj_val * value
                if prod == -1:
                    ndifferent_next += 1
                elif prod == 1:
                    nsame_next += 1
                else:
                    raise ValueError('edge sign error')
        # If the target number of edges with sign changes
        # is compatible with the current known edge change status
        # then we are OK.
        if is_branch_compat(nsame_next, ndifferent_next, ntarget, nbranches):
            id_to_val[idcur] = value
            if depth == len(internals) - 1:
                yield dict(id_to_val)
            else:
                for v in rec_internal(
                        id_to_adj, id_to_val,
                        nsame_next, ndifferent_next, ntarget, nbranches,
                        internals, depth+1):
                    yield v
        # Reset the current value to None.
        id_to_val[idcur] = None

def gen_assignments(
        id_to_adj, id_to_val,
        ntarget, nbranches, internals):
    """
    This is the facade for a recursive function.
    """
    # define the parameters for the recursive function
    # create the generator object
    nsame = 0
    ndifferent = 0
    depth = 0
    obj = rec_internal(
            id_to_adj, id_to_val,
            nsame, ndifferent, ntarget, nbranches,
            internals, depth)
    # return the generator object
    return obj

def gen_assignments_brute(
        id_to_adj, id_to_val,
        ntarget, nbranches, internals):
    """
    This is like gen_assignments except it is slow.
    This is for testing.
    """
    ninternals = len(internals)
    for values in itertools.product((-1, 1), repeat=ninternals):
        for x, value in zip(internals, values):
            id_to_val[x] = value
        id_to_region = get_regions(id_to_adj, id_to_val)
        nregions = len(set(id_to_region.values()))
        if nregions == ntarget + 1:
            yield dict(id_to_val)
    for x in internals:
        id_to_val[x] = None

def rec_eigen_weak(id_to_adj, id_to_val_list, id_to_list_val,
        depth, require_sign_harmonicity):
    """
    This is a recursive function.
    Each level corresponds to an eigenvector.
    This uses the relatively weak condition of principal orthant connectivity.
    @param id_to_adj: maps an id to a list of adjacent ids
    @param id_to_val_list: a list of k partial valuation maps
    @param id_to_list_val: maps an id to a list of values
    @param depth: zero corresponds to fiedler depth
    @param require_sign_harmonicity: True if sign harmonicity is required
    @return: None or a valid map
    """
    # Define the set of ids.
    ids = set(id_to_adj)
    # Define the requested number of cut branches at this depth.
    ntarget = depth + 1
    # Get the number of branches in the tree.
    nbranches = sum(len(v) for v in id_to_adj.values()) / 2
    # Get the list of internal ids.
    internals = sorted(get_internal_set(id_to_adj))
    # Consider each assignment at this level that satisfies ntarget.
    for d in gen_assignments(
            id_to_adj, id_to_val_list[depth],
            ntarget, nbranches, internals):
        # Check sign harmonicity if requested.
        if require_sign_harmonicity:
            if not is_sign_harmonic(id_to_adj, d):
                continue
        # make the putative next cumulative valuation
        id_to_list_next = {}
        for x in ids:
            v = id_to_list_val.get(x, [])
            id_to_list_next[x] = tuple(list(v) + [d[x]])
        # Require the cumulative assignment to meet orthant connectivity.
        if is_value_connected(id_to_adj, id_to_list_next):
            if depth == len(id_to_val_list) - 1:
                return id_to_list_next
            else:
                result = rec_eigen_weak(
                        id_to_adj, id_to_val_list, id_to_list_next,
                        depth + 1, require_sign_harmonicity)
                if result:
                    return result

def rec_eigen_strong(id_to_adj, id_to_val_list, id_to_list_val,
        depth, require_sign_harmonicity):
    """
    This is a recursive function.
    Each level corresponds to an eigenvector.
    This uses the stronger condition relating sign graphs.
    @param id_to_adj: maps an id to a list of adjacent ids
    @param id_to_val_list: a list of k partial valuation maps
    @param id_to_list_val: maps an id to a list of values
    @param depth: zero corresponds to fiedler depth
    @param require_sign_harmonicity: True if sign harmonicity is required
    @return: None or a valid map
    """
    # Define the set of ids.
    ids = set(id_to_adj)
    # Define the requested number of cut branches at this depth.
    ntarget = depth + 1
    # Get the number of branches in the tree.
    nbranches = sum(len(v) for v in id_to_adj.values()) / 2
    # Get the list of internal ids.
    internals = sorted(get_internal_set(id_to_adj))
    # Consider each assignment at this level that satisfies ntarget.
    for d in gen_assignments(
            id_to_adj, id_to_val_list[depth],
            ntarget, nbranches, internals):
        # Check sign harmonicity if requested.
        if require_sign_harmonicity:
            if not is_sign_harmonic(id_to_adj, d):
                continue
        # make the putative next cumulative valuation
        id_to_list_next = {}
        for x in ids:
            v = id_to_list_val.get(x, [])
            id_to_list_next[x] = tuple(list(v) + [d[x]])
        # get the valuation at the previous depth
        if depth:
            d_prev = dict((x, id_to_list_val[x][-1]) for x in ids)
        else:
            d_prev = dict((x, 1) for x in ids)
        # Check the pairwise sign lacing condition.
        if check_sign_lacing(id_to_adj, d_prev, d):
            if depth == len(id_to_val_list) - 1:
                return id_to_list_next
            else:
                result = rec_eigen_strong(
                        id_to_adj, id_to_val_list, id_to_list_next,
                        depth + 1, require_sign_harmonicity)
                if result:
                    return result

def check_sign_lacing(id_to_adj, id_to_va, id_to_vb):
    """
    @param id_to_adj: maps an id to a list of adjacent ids
    @param id_to_va: map id to sign valuation k
    @param id_to_vb: map id to sign valuation k+1
    @return: True if the interlacing is compatible
    """
    ids = set(id_to_adj)
    id_to_region = get_regions(id_to_adj, id_to_va)
    id_to_pair = dict((x, (id_to_region[x], id_to_vb[x])) for x in ids)
    return is_value_connected(id_to_adj, id_to_pair)

def get_internal_set(id_to_adj):
    return set(v for v, d in id_to_adj.items() if len(d) > 1)

def get_leaf_set(id_to_adj):
    return set(v for v, d in id_to_adj.items() if len(d) == 1)

def get_leaf_lists(id_to_adj, id_to_val):
    """
    Find leaf ids with common values.
    @param id_to_adj: maps an id to a list of adjacent ids
    @param id_to_val: maps an id to a value
    @return: a list of lists of leaf ids
    """
    value_to_set = defaultdict(set)
    leaves = get_leaf_set(id_to_adj)
    for leaf in leaves:
        val = id_to_val[leaf]
        value_to_set[val].add(leaf)
    return [list(s) for s in value_to_set.values()]

def is_value_connected(id_to_adj, id_to_val):
    """
    Note that in this case the value can be a tuple of values.
    This function checks that in the tree,
    elements with the same values are connected.
    Note that id_to_adj is assumed to be a tree,
    and the values in id_to_val must be hashable.
    Here are two usage examples.
    The first usage example is to check principal orthant connectivity
    by looking at the tuple of the first k valuations.
    The second usage example is to check a stronger
    connectivity criterion by looking at the pair
    such that the first element is a region index for the kth valuation
    and where the second element is the (k+1)st valuation itself.
    @param id_to_adj: maps an id to a list of adjacent ids
    @param id_to_val: maps an id to a value
    """
    id_to_region = get_regions(id_to_adj, id_to_val)
    nvalues = len(set(id_to_val.values()))
    nregions = len(set(id_to_region.values()))
    return nvalues == nregions

def get_regions(id_to_adj, id_to_val):
    """
    Find connected regions with uniform value.
    Assume a tree topology.
    Each region will get an arbitrary color.
    @param id_to_adj: maps an id to a list of adjacent ids
    @param id_to_val: maps an id to a value
    @return: a map from id to region
    """
    # begin with the min id for determinism for testing
    x = min(id_to_adj)
    id_to_region = {x : 0}
    nregions = 1
    shell = set([x])
    visited = set([x])
    while shell:
        next_shell = set()
        for v in shell:
            v_val = id_to_val[v]
            v_region = id_to_region[v]
            # sort for determinism for testing
            for u in sorted(id_to_adj[v]):
                if u not in visited:
                    u_val = id_to_val[u]
                    if u_val == v_val:
                        id_to_region[u] = v_region
                    else:
                        id_to_region[u] = nregions
                        nregions += 1
                    visited.add(u)
                    next_shell.add(u)
        shell = next_shell
    return id_to_region

def is_sign_harmonic(id_to_adj, id_to_val):
    """
    Sign harmonic will mean that each strong sign graph has a leaf.
    Assume all values are either +1 or -1.
    @param id_to_adj: maps an id to a list of adjacent ids
    @param id_to_val: maps an id to a value
    """
    leaves = get_leaf_set(id_to_adj)
    visited = set(leaves)
    shell = set(leaves)
    while shell:
        next_shell = set()
        for v in shell:
            v_val = id_to_val[v]
            for u in id_to_adj[v]:
                if u not in visited:
                    u_val = id_to_val[u]
                    if u_val == v_val:
                        visited.add(u)
                        next_shell.add(u)
        shell = next_shell
    nvertices = len(id_to_adj)
    nvisited = len(visited)
    return nvertices == nvisited


g_test_id_to_adj = {
        1 : [6],
        2 : [6],
        3 : [8],
        4 : [7],
        5 : [7],
        6 : [1, 2, 8],
        7 : [4, 5, 8],
        8 : [3, 6, 7]}

def d_to_fset_of_pairs(d):
    """
    @param d: a dict
    @return: a frozenset of pairs
    """
    return frozenset((k, v) for k, v in d.items())

class TestThis(unittest.TestCase):

    def test_gen_assignments_a(self):
        id_to_val = {1: 1, 2: 1, 3: 1, 4: 1, 5: -1, 6: None, 7: None, 8: None}
        ntarget = 2
        nbranches = 7
        internals = [6, 7, 8]
        # Get the list of observed compatible assignments.
        ds = list(gen_assignments(
                g_test_id_to_adj, id_to_val,
                ntarget, nbranches, internals))
        # Define the only true compatible assignment.
        d = {1: 1, 2: 1, 3: 1, 4: 1, 5: -1, 6: 1, 7: -1, 8: 1}
        # Compare the observed and expected assignments.
        self.assertEqual(len(ds), 1)
        self.assertEqual(ds[0], d)

    def test_gen_assignments_b(self):
        id_to_adj = {
                1: [5],
                2: [5],
                3: [6],
                4: [6],
                5: [1, 2, 6],
                6: [3, 4, 5]}
        id_to_val_list = [
                {1:-1, 2:-1, 3:1, 4:1, 5:None, 6:None},
                {1:-1, 2:1, 3:1, 4:-1, 5:None, 6:None},
                {1:-1, 2:-1, 3:1, 4:-1, 5:None, 6:None}]
        for i, id_to_val in enumerate(id_to_val_list):
            ntarget = i + 1
            nbranches = 5
            internals = [5, 6]
            fast_list_of_dicts = list(gen_assignments(
                id_to_adj, id_to_val,
                ntarget, nbranches, internals))
            slow_list_of_dicts = list(gen_assignments_brute(
                id_to_adj, id_to_val,
                ntarget, nbranches, internals))
            fast_frozensets = set(
                    d_to_fset_of_pairs(d) for d in fast_list_of_dicts)
            slow_frozensets = set(
                    d_to_fset_of_pairs(d) for d in slow_list_of_dicts)
            # The sets should not be empty.
            self.assertTrue(len(fast_frozensets) > 0)
            self.assertTrue(len(slow_frozensets) > 0)
            # The sets should be equal.
            self.assertEqual(fast_frozensets, slow_frozensets)

    def test_sign_harmonic_a(self):
        id_to_val = {
                1:1, 2:1, 3:1, 4:1, 5:1, 6:1, 7:1, 8:1}
        observed = is_sign_harmonic(g_test_id_to_adj, id_to_val)
        expected = True
        self.assertEqual(observed, expected)

    def test_sign_harmonic_b(self):
        id_to_val = {
                1:-1, 2:-1, 3:1, 4:1, 5:1, 6:-1, 7:1, 8:1}
        observed = is_sign_harmonic(g_test_id_to_adj, id_to_val)
        expected = True
        self.assertEqual(observed, expected)

    def test_sign_harmonic_c(self):
        id_to_val = {
                1:1, 2:1, 3:1, 4:1, 5:1, 6:-1, 7:1, 8:-1}
        observed = is_sign_harmonic(g_test_id_to_adj, id_to_val)
        expected = False
        self.assertEqual(observed, expected)

    def test_get_regions_a(self):
        id_to_val = {
                1:1, 2:1, 3:1, 4:1, 5:1, 6:-1, 7:1, 8:-1}
        observed = get_regions(g_test_id_to_adj, id_to_val)
        expected = {1:0, 6:1, 2:2, 8:1, 3:3, 7:4, 5:4, 4:4}
        self.assertEqual(observed, expected)

    def test_get_regions_b(self):
        id_to_val = {
                1:1, 2:1, 3:1, 4:1, 5:1, 6:1, 7:1, 8:1}
        observed = get_regions(g_test_id_to_adj, id_to_val)
        expected = {1:0, 2:0, 3:0, 4:0, 5:0, 6:0, 7:0, 8:0}
        self.assertEqual(observed, expected)

    def test_orthant_connected_a(self):
        id_to_val = {
                1 : (1, 1),
                2 : (1, 1),
                3 : (-1, -1),
                4 : (-1, -1),
                5 : (-1, 1),
                6 : (1, 1),
                7 : (-1, -1),
                8 : (-1, -1)}
        observed = is_value_connected(g_test_id_to_adj, id_to_val)
        expected = True
        self.assertEqual(observed, expected)

    def test_orthant_connected_b(self):
        id_to_val = {
                1 : (1, 1),
                2 : (1, 1),
                3 : (-1, -1),
                4 : (-1, -1),
                5 : (-1, 1),
                6 : (1, 1),
                7 : (-1, 1),
                8 : (-1, -1)}
        observed = is_value_connected(g_test_id_to_adj, id_to_val)
        expected = False
        self.assertEqual(observed, expected)

    def test_rec_eigen_weak_a(self):
        id_to_val_list = [
                {1:1, 2:1, 3:-1, 4:-1, 5:-1, 6:None, 7:None, 8:None},
                {1:1, 2:1, 3:-1, 4:1, 5:1, 6:None, 7:None, 8:None}]
        id_to_list_val = {}
        observed = rec_eigen_weak(
                g_test_id_to_adj, id_to_val_list, id_to_list_val,
                0, True)
        expected = {
                1: (1, 1),
                2: (1, 1),
                3: (-1, -1),
                4: (-1, 1),
                5: (-1, 1),
                6: (1, 1),
                7: (-1, 1),
                8: (-1, -1)}
        self.assertEqual(observed, expected)

    def test_rec_eigen_weak_b(self):
        id_to_val_list = [
                {1:1, 2:1, 3:-1, 4:-1, 5:-1, 6:None, 7:None, 8:None},
                {1:1, 2:1, 3:1, 4:1, 5:-1, 6:None, 7:None, 8:None}]
        id_to_list_val = {}
        observed = rec_eigen_weak(
                g_test_id_to_adj, id_to_val_list, id_to_list_val,
                0, True)
        expected = None
        self.assertEqual(observed, expected)

    def test_rec_eigen_strong(self):
        id_to_adj = {
                1: [5],
                2: [5],
                3: [6],
                4: [6],
                5: [1, 2, 6],
                6: [3, 4, 5]}
        id_to_val_list = [
                {1:-1, 2:-1, 3:1, 4:1, 5:None, 6:None},
                {1:-1, 2:1, 3:1, 4:-1, 5:None, 6:None},
                {1:-1, 2:-1, 3:1, 4:-1, 5:None, 6:None}]
        id_to_list_val = {}
        observed = rec_eigen_strong(
                id_to_adj, id_to_val_list, id_to_list_val,
                0, True)
        expected = {
                1: (-1, -1, -1),
                2: (-1, 1, -1),
                3: (1, 1, 1),
                4: (1, -1, -1),
                5: (-1, 1, 1),
                6: (1, 1, 1)}
        self.assertEqual(observed, expected)

    def test_check_sign_lacing_true(self):
        id_to_adj = {
                1: [5],
                2: [5],
                3: [6],
                4: [6],
                5: [1, 2, 6],
                6: [3, 4, 5]}
        vs = [
                {1:1, 2:1, 3:1, 4:1, 5:1, 6:1},
                {1:-1, 2:-1, 3:1, 4:1, 5:-1, 6:1},
                {1:-1, 2:1, 3:1, 4:-1, 5:1, 6:1},
                {1:-1, 2:-1, 3:1, 4:-1, 5:1, 6:1}]
        for va, vb in iterutils.pairwise(vs):
            observed = check_sign_lacing(id_to_adj, va, vb)
            expected = True
            self.assertEqual(observed, expected)

    def test_check_sign_lacing_false(self):
        id_to_adj = {
                1: [5],
                2: [5],
                3: [6],
                4: [6],
                5: [1, 2, 6],
                6: [3, 4, 5]}
        va = {1:-1, 2:1, 3:1, 4:-1, 5:1, 6:1}
        vb = {1:-1, 2:1, 3:1, 4:-1, 5:1, 6:-1}
        observed = check_sign_lacing(id_to_adj, va, vb)
        expected = False
        self.assertEqual(observed, expected)


if __name__ == '__main__':
    unittest.main()

