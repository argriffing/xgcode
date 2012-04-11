"""
Split algebra module.

This defines some split algebra used by Ben Redelings
to study what he calls multiconnected trees.
Most things in this module are frozensets of frozensets or integers.
The things referred to as states are integer indices starting from zero.
In a phylogenetic setting the states are leaf taxon indices.
"""

import unittest
import itertools
from itertools import combinations
from itertools import permutations

def mksplit(a, b):
    return frozenset((frozenset(a), frozenset(b)))

def is_informative(split):
    return min(len(side) for side in split) > 1

def get_quartets(N):
    if N < 4:
        return set()
    states = set(range(N))
    accum = set()
    for i, j, k, l in combinations(states, 4):
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
    states = set(range(N))
    accum = set()
    for a_size in range(2, N - 1):
        for a in combinations(states, a_size):
            accum.add(mksplit(a, states - set(a)))
    return accum

def get_bifurcating_trees(N):
    """
    A bifurcating tree is a frozenset of jointly compatible full splits.
    It has exactly N-3 such splits.
    """
    if N < 3:
        raise ValueError('not enough states')
    states = set(range(N))
    trees = set()
    for splits in combinations(get_informative_full_splits(N), N-3):
        if pairwise_split_compatibility(splits):
            trees.add(frozenset(splits))
    return trees

def get_multifurcating_trees(N):
    """
    A bifurcating tree is a frozenset of jointly compatible full splits.
    It has at most N-3 such splits.
    """
    if N < 3:
        raise ValueError('not enough states')
    states = set(range(N))
    all_splits = get_informative_full_splits(N)
    trees = set()
    for nsplits in range(N-2):
        for splits in combinations(all_splits, nsplits):
            if pairwise_split_compatibility(splits):
                trees.add(frozenset(splits))
    return trees

def get_mc_trees_slow(N):
    """
    A multiconnected tree is a frozenset of jointly compatible partial splits.
    It has at most N-3 such splits.
    It is also subject to other constraints.
    """
    if N < 3:
        raise ValueError('not enough states')
    states = set(range(N))
    all_splits = get_informative_partial_splits(N)
    trees = set()
    for nsplits in range(N-2):
        for splits in combinations(all_splits, nsplits):
            if pairwise_split_consistency(splits):
                trees.add(frozenset(splits))
    return trees

def get_mc_trees(N):
    """
    A multiconnected tree is a frozenset of jointly compatible partial splits.
    It has at most N-3 such splits.
    It is also subject to other constraints, namely consistency.
    """
    if N < 3:
        raise ValueError('not enough states')
    states = set(range(N))
    all_splits = list(get_informative_partial_splits(N))
    nsplits_total = len(all_splits)
    # precompute the pairwise consistencies
    d = {}
    for i, split_i in enumerate(all_splits):
        for j, split_j in enumerate(all_splits):
            if i < j:
                d[i, j] = split_consistency(split_i, split_j)
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
            tree = frozenset(all_splits[i] for i in index_set)
            accum_trees.add(tree)
    return accum_trees

def get_nearly_mc_trees(N):
    """
    This tree type is also a frozenset of jointly compatible partial splits.
    It has at most N-3 such splits.
    Adds the criteria of compatibility and non-redundancy but not consistency.
    """
    if N < 3:
        raise ValueError('not enough states')
    states = set(range(N))
    all_splits = list(get_informative_partial_splits(N))
    nsplits_total = len(all_splits)
    # precompute the pairwise agreeability
    d = {}
    for i, a in enumerate(all_splits):
        for j, b in enumerate(all_splits):
            if i < j:
                d[i, j] = all([
                    split_compatibility(a, b),
                    not split_redundancy(a, b)])
    # --- the following is copypasted from mc-trees ---
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
            tree = frozenset(all_splits[i] for i in index_set)
            accum_trees.add(tree)
    return accum_trees

def split_implication(split_a, split_b):
    """
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

def split_consistency(split_a, split_b, verbose=False):
    result, msg = _split_consistency(split_a, split_b)
    if verbose:
        return result, msg
    else:
        return result

def pairwise_split_compatibility(splits):
    return all(split_compatibility(a,b) for a,b in combinations(splits, 2))

def pairwise_split_redundancy(splits):
    return any(split_redundancy(a,b) for a,b in combinations(splits, 2))

def pairwise_split_consistency(splits):
    return all(split_consistency(a,b) for a,b in combinations(splits, 2))

def tree_compatibility(tree_a, tree_b):
    """
    For the purpose of this function, a tree is an iterable of splits.
    More stringently, a tree should be a frozenset of jointly compatible splits.
    """
    return all(split_compatibility(a, b) for a in tree_a for b in tree_b)

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

    def test_get_multifurcating_trees(self):
        """
        Sloane sequence A000311.
        """
        expected = [4, 26, 236, 2752]
        Ns = range(4, 4 + len(expected))
        observed = [len(get_multifurcating_trees(N)) for N in Ns]
        self.assertEqual(observed, expected)

    def test_get_mc_trees(self):
        """
        Sloane sequence does not exist.
        The extended sequence is [4, 41, 746, 20462, 780886, ...].
        """
        expected = [4, 41, 746]
        Ns = range(4, 4 + len(expected))
        observed = [len(get_mc_trees(N)) for N in Ns]
        self.assertEqual(observed, expected)

    def test_get_nearly_mc_trees(self):
        """
        Sloane sequence does not exist.
        The extended sequence is [4, 191, 124186, ...].
        """
        expected = [4, 1199]
        Ns = range(4, 4 + len(expected))
        observed = [len(get_nearly_mc_trees(N)) for N in Ns]
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
