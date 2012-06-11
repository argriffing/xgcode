"""
Split algebra module.

This defines some split algebra used by Ben Redelings
to study what he calls multiconnected trees.
Most things in this module are frozensets of frozensets or integers.
"""

import unittest
import itertools
from itertools import combinations
from itertools import permutations

import Util
import iterutils

def mksplit(a, b):
    return frozenset((frozenset(a), frozenset(b)))

def is_informative(split):
    return min(len(side) for side in split) > 1

def get_quartets(N):
    if N < 4:
        return set()
    accum = set()
    for i, j, k, l in combinations(range(N), 4):
        accum.add(mksplit((i,j), (k,l)))
        accum.add(mksplit((i,k), (j,l)))
        accum.add(mksplit((i,l), (j,k)))
    return accum

def get_informative_partial_splits(N):
    if N < 4:
        return set()
    accum = set()
    for remaining_size in range(4, N+1):
        for remaining in combinations(range(N), remaining_size):
            remaining_set = set(remaining)
            for a_size in range(2, remaining_size - 1):
                for a in combinations(remaining, a_size):
                    accum.add(mksplit(a, remaining_set - set(a)))
    return accum

def get_informative_full_splits(N):
    if N < 4:
        return set()
    accum = set()
    for a_size in range(2, N - 1):
        for a in combinations(range(N), a_size):
            accum.add(mksplit(a, set(range(N)) - set(a)))
    return accum

def get_bifurcating_trees(N):
    """
    A bifurcating tree is a frozenset of N-3 pairwise compatible full splits.
    """
    trees = set()
    for splits in combinations(get_informative_full_splits(N), N-3):
        if pairwise_split_compatibility(splits):
            trees.add(frozenset(splits))
    return trees

def count_bifurcating_trees(N):
    p = 1
    for i in range(1, N-2):
        p *= 2*i + 1
    return p

def get_multifurcating_trees(N):
    """
    A multifurcating tree is a frozenset of pairwise compatible full splits.
    """
    all_splits = get_informative_full_splits(N)
    trees = set()
    for nsplits in range(N-2):
        shell = set()
        for splits in combinations(all_splits, nsplits):
            if pairwise_split_compatibility(splits):
                shell.add(frozenset(splits))
        if not shell:
            break
        trees.update(shell)
    return trees

def count_multifurcating_trees(N):
    a = [0, 1, 1]
    for n in range(2, N-1):
        sigma = sum(Util.choose(n, k)*a[k]*a[n-k+1] for k in range(2, n))
        b = (n+2) * a[n] + 2*sigma
        a.append(b)
    return a[-1]

def count_independently_unplugged_trees(N):
    c = 1
    for n_known in range(4, N+1):
        n_at_large = N - n_known
        coeff = Util.choose(N, n_at_large)
        c += coeff * (count_multifurcating_trees(n_known) - 1)
    return c

def get_mc_trees_slow(N):
    """
    A multiconnected tree is a frozenset of pairwise compatible partial splits.
    It is also subject to other constraints, namely consistency.
    """
    all_splits = get_informative_partial_splits(N)
    trees = set()
    for nsplits in itertools.count():
        shell = set()
        for splits in combinations(all_splits, nsplits):
            if pairwise_split_consistency(splits):
                shell.add(frozenset(splits))
        if not shell:
            break
        trees.update(shell)
    return trees

def get_mc_trees(N):
    """
    A multiconnected tree is a frozenset of pairwise compatible partial splits.
    It is also subject to other constraints, namely consistency.
    """
    partial_splits = list(get_informative_partial_splits(N))
    d = {}
    for i, a in enumerate(partial_splits):
        for j, b in enumerate(partial_splits):
            if i < j:
                d[i, j] = split_consistency(a, b)
    return _get_pairwise_admissible_trees(d, partial_splits)

def get_quartet_sets(N):
    """
    The quartets within each set are pairwise compatible.
    This should be a subset of nearly_mc_trees.
    """
    quartets = list(get_quartets(N))
    d = {}
    for i, a in enumerate(quartets):
        for j, b in enumerate(quartets):
            if i < j:
                d[i, j] = split_compatibility(a, b)
    return _get_pairwise_admissible_trees(d, quartets)

def get_nearly_mc_trees(N):
    """
    This tree type is also a frozenset of pairwise compatible partial splits.
    Adds the criteria of compatibility and non-redundancy but not consistency.
    """
    partial_splits = list(get_informative_partial_splits(N))
    d = {}
    for i, a in enumerate(partial_splits):
        for j, b in enumerate(partial_splits):
            if i < j:
                d[i, j] = all([
                    split_compatibility(a, b),
                    not split_redundancy(a, b)])
    return _get_pairwise_admissible_trees(d, partial_splits)

def _get_pairwise_admissible_trees(d, partial_splits):
    """
    Get all split combinations that meet the pairwise admissibility criterion.
    @param d: pairwise index admissibility dictionary
    @param partial_splits: list of all informative partial splits
    @return: a set of trees where each tree is a frozenset of partial splits
    """
    nsplits_total = len(partial_splits)
    accum_trees = set([frozenset([])])
    next_shell = set([frozenset([])])
    while next_shell:
        shell = next_shell
        next_shell = set()
        for index_set in shell:
            for j in range(1+max(index_set | frozenset([-1])), nsplits_total):
                if all(d[i, j] for i in index_set):
                    next_index_set = index_set | frozenset([j])
                    next_shell.add(next_index_set)
        for index_set in next_shell:
            tree = frozenset(partial_splits[i] for i in index_set)
            accum_trees.add(tree)
    return accum_trees

def get_representable_treesets(N, pre_trees):
    """
    Get a set of representations of p-representable fully resolved tree sets.
    A p-representable set is a set of fully resolved trees
    that is representable by a set of partial splits.
    The pre_trees parameter could for example be the set of
    pairwise compatible and pairwise-non-redundant partial splits.
    Or it could be the set of pairwise consistent partial splits.
    @param N: the number of tips
    @param pre_trees: a candidate sequence of frozensets of partial splits
    @return: set of sets of bifurcating trees
    """
    treeset_sets = set()
    for pre_tree in pre_trees:
        trees = list(get_bifurcating_trees(N))
        # repeatedly cut down the set of compatible trees
        for split in pre_tree:
            trees = [t for t in trees if split_tree_compatibility(split, t)]
        # if a set of trees is compatible then add it
        if trees:
            treeset_sets.add(frozenset(trees))
    return treeset_sets

def resolved_to_quartets(resolved):
    """
    @param resolved: a bifurcating tree defined by compatible full splits
    @return: set of compatible quartets
    """
    quartets = set()
    for split in resolved:
        a, b = split
        for qa in combinations(a, 2):
            for qb in combinations(b, 2):
                quartets.add(mksplit(qa, qb))
    return quartets

def split_implication(split_a, split_b):
    """
    Defined by Ben Redelings.
    This split relation is asymmetric.
    For example, 012|34 implies 12|34.
    Also, 12|34 implies 12|34.
    @param split_a: a partial split
    @param split_b: a partial split
    @return: True if split_a implies split_b
    """
    a1, a2 = split_a
    b1, b2 = split_b
    return (b1 <= a1 and b2 <= a2) or (b2 <= a1 and b1 <= a2)

def split_redundancy(a, b):
    """
    This split relation is symmetric.
    """
    return split_implication(a, b) or split_implication(b, a)

def split_compatibility(split_a, split_b):
    """
    This split relation is symmetric.
    """
    a1, a2 = split_a
    b1, b2 = split_b
    return not all((a1 & b1, a1 & b2, a2 & b1, a2 & b2))

def _split_consistency(split_a, split_b):
    """
    This split relation is symmetric.
    True if the pair is consistent with at least one of the nine conditions.
    """
    a1, a2 = split_a
    b1, b2 = split_b
    at2, at1 = split_a
    bt2, bt1 = split_b
    A = (a1, a2)
    AT = (at1, at2)
    B = (b1, b2)
    BT = (bt1, bt2)
    # check the 'left of' conditions
    if a1 < b1 and b2 < a2:
        return True, '%s < %s' % (A, B)
    if b1 < a1 and a2 < b2:
        return True, '%s < %s' % (B, A)
    if at1 < b1 and b2 < at2:
        return True, '%s < %s' % (AT, B)
    if b1 < at1 and at2 < b2:
        return True, '%s < %s' % (B, AT)
    # check the 'wanders over' conditions
    if b1 | b2 <= a2:
        return True, '%s wanders over %s' % (A, B)
    if b1 | b2 <= at2:
        return True, '%s wanders over %s' % (AT, B)
    if a1 | a2 <= b2:
        return True, '%s wanders over %s' % (B, A)
    if a1 | a2 <= bt2:
        return True, '%s wanders over %s' % (BT, A)
    # check the orthogonality condition
    if not ((a1 | a2) & (b1 | b2)):
        return True, '%s is orthogonal to %s' % (A, B)
    return False, 'fails the rules'

def split_consistency_verbose(split_a, split_b):
    return _split_consistency(split_a, split_b)

def split_consistency(split_a, split_b):
    result, message = _split_consistency(split_a, split_b)
    return result

def pairwise_split_compatibility(splits):
    return all(split_compatibility(a,b) for a,b in combinations(splits, 2))

def pairwise_split_redundancy(splits):
    return any(split_redundancy(a,b) for a,b in combinations(splits, 2))

def pairwise_split_consistency(splits):
    return all(split_consistency(a,b) for a,b in combinations(splits, 2))

def split_tree_compatibility(split, tree):
    return all(split_compatibility(split, s) for s in tree)

def all_quartets_are_individually_implied(implied_quartets, partial_splits):
    """
    This function checks the following statement.
    For all implied quartets there is a partial split that implies that quartet.
    @param implied_quartets: quartets implied by a tree
    @param partial_splits: partial splits
    """
    qs = implied_quartets
    ps = partial_splits
    return all(any(split_implication(p, q) for p in ps) for q in qs)




class TestSplitAlgebra(unittest.TestCase):

    def test_get_quartets(self):
        """
        Sloane sequence A050534.
        """
        expected = [3, 15, 45, 105, 210, 378, 630, 990, 1485]
        Ns = range(4, 4 + len(expected))
        observed = [len(get_quartets(N)) for N in Ns]
        self.assertEqual(observed, expected)

    def test_get_bifurcating_trees(self):
        """
        Sloane sequence A001147.
        """
        expected = [3, 15, 105, 945]
        Ns = range(4, 4 + len(expected))
        observed = [len(get_bifurcating_trees(N)) for N in Ns]
        self.assertEqual(observed, expected)

    def test_count_bifurcating_trees(self):
        """
        Sloane sequence A001147.
        """
        expected = [3, 15, 105, 945, 10395, 135135, 2027025, 34459425]
        Ns = range(4, 4 + len(expected))
        observed = [count_bifurcating_trees(N) for N in Ns]
        self.assertEqual(observed, expected)

    def test_get_multifurcating_trees(self):
        """
        Sloane sequence A000311.
        """
        expected = [4, 26, 236, 2752]
        Ns = range(4, 4 + len(expected))
        observed = [len(get_multifurcating_trees(N)) for N in Ns]
        self.assertEqual(observed, expected)

    def test_count_multifurcating_trees(self):
        """
        Sloane sequence A000311.
        """
        expected = [4, 26, 236, 2752]
        Ns = range(4, 4 + len(expected))
        observed = [count_multifurcating_trees(N) for N in Ns]
        self.assertEqual(observed, expected)

    def test_count_independently_unplugged_trees(self):
        """
        ???
        """
        expected = [4, 41, 431, 5027, 69406, 1135199, 21569937]
        Ns = range(4, 4 + len(expected))
        observed = [count_independently_unplugged_trees(N) for N in Ns]
        self.assertEqual(observed, expected)

    def test_get_mc_trees_slow(self):
        """
        Sloane sequence does not exist.
        The extended sequence is [4, 41, 746, 20462, 780886, ???].
        """
        expected = [4, 41]
        Ns = range(4, 4 + len(expected))
        observed = []
        for N in Ns:
            trees = get_mc_trees_slow(N)
            self.assertIn(frozenset([]), trees)
            observed.append(len(trees))
        self.assertEqual(observed, expected)

    def test_get_mc_trees(self):
        """
        Sloane sequence does not exist.
        The extended sequence is [4, 41, 746, 20462, 780886, ???].
        """
        expected = [4, 41, 746]
        Ns = range(4, 4 + len(expected))
        observed = []
        for N in Ns:
            trees = get_mc_trees(N)
            self.assertIn(frozenset([]), trees)
            observed.append(len(trees))
        self.assertEqual(observed, expected)

    def test_get_nearly_mc_trees(self):
        """
        Sloane sequence does not exist.
        The extended sequence is [4, 1199, ???].
        """
        expected = [4, 1199]
        Ns = range(4, 4 + len(expected))
        observed = [len(get_nearly_mc_trees(N)) for N in Ns]
        self.assertEqual(observed, expected)

    def test_get_quartet_sets(self):
        """
        The extended sequence is [4, 1024, ???].
        """
        expected = [4, 1024]
        Ns = range(4, 4 + len(expected))
        observed = [len(get_quartet_sets(N)) for N in Ns]
        self.assertEqual(observed, expected)

    def test_get_representable_treesets_p(self):
        """
        The extended sequence is [4, 41, ???].
        """
        expected = [4, 41]
        Ns = range(4, 4 + len(expected))
        observed = []
        for N in Ns:
            pre_trees = list(get_nearly_mc_trees(N))
            observed.append(len(get_representable_treesets(N, pre_trees)))
        self.assertEqual(observed, expected)

    def test_get_representable_treesets_q(self):
        """
        The extended sequence is [4, 41, ???].
        """
        expected = [4, 41]
        Ns = range(4, 4 + len(expected))
        observed = []
        for N in Ns:
            pre_trees = list(get_quartet_sets(N))
            observed.append(len(get_representable_treesets(N, pre_trees)))
        self.assertEqual(observed, expected)

    def test_quartet_set_representable_tree_sets(self):
        """
        Sloane sequence does not exist.
        The extended sequence is [4, 41, 1586, ???].
        """
        expected = [4, 41]
        Ns = range(4, 4 + len(expected))
        observed = []
        for N in Ns:
            # get all quartets
            qs = list(get_quartets(N))
            # get all resolved trees
            ts = list(get_bifurcating_trees(N))
            # map from quartet to quartet index
            q_to_i = dict((q, i) for i, q in enumerate(qs))
            # map tree index to set of indices of compatible quartets
            ti_to_qi_set = dict((i, set(q_to_i[q] for q in resolved_to_quartets(t))) for i, t in enumerate(ts))
            # begin the definition of representable tree index sets
            representable_ti_sets = set()
            # For each tree look at every subset of compatible quartets.
            # Record the set of tree indices compatible with each subset.
            for t in ts:
                tqis = [q_to_i[q] for q in resolved_to_quartets(t)]
                for qi_subset_tuple in iterutils.powerset(tqis):
                    qi_subset = set(qi_subset_tuple)
                    #ti_pattern = frozenset(i for i, tree in enumerate(ts) if 
                    ti_pattern = frozenset(i for i, s in ti_to_qi_set.items() if qi_subset <= s)
                    representable_ti_sets.add(ti_pattern)
            n = len(representable_ti_sets)
            observed.append(n)
        self.assertEqual(observed, expected)

    def test_get_informative_partial_splits(self):
        """
        Sloane sequence A112495.
        """
        expected = [3, 25, 130, 546, 2037]
        Ns = range(4, 4 + len(expected))
        observed = [len(get_informative_partial_splits(N)) for N in Ns]
        self.assertEquals(observed, expected)

    def test_get_informative_full_splits(self):
        """
        Sloane sequence A000247.
        """
        expected = [3, 10, 25, 56, 119]
        Ns = range(4, 4 + len(expected))
        observed = [len(get_informative_full_splits(N)) for N in Ns]
        self.assertEquals(observed, expected)

    def test_ben_counterexample_a(self):
        """
        This counterexample is from a manuscript.
        """
        a = mksplit([0,1], [2,3])
        b = mksplit([0,1], [4,5])
        self.assertFalse(split_consistency(a, b))

    def test_ben_counterexample_b(self):
        """
        This counterexample is from a manuscript.
        """
        a = mksplit([0,1,6], [2,3,4,5])
        b = mksplit([0,1,2,3], [4,5])
        self.assertFalse(split_consistency(a, b))

    def test_ben_example_a(self):
        """
        This counter-counterexample is from a manuscript.
        """
        a = mksplit([0,1,6], [2,3,4,5])
        b = mksplit([0,1,2,3,6], [4,5])
        self.assertTrue(split_consistency(a, b))


if __name__ == '__main__':
    unittest.main()
