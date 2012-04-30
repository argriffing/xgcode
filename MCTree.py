"""
Multiconnected tree module.

This is a concept by Ben Redelings in
an unpublished manuscript (as of June 2009 - April 2012).
This module aims to explore some of the combinatorics
of multiconnected trees.
"""

import unittest

import splitalg
from splitalg import mksplit


# TODO move the content of the "tests" into the command line of web scripts
# because they are not actually tests and they are slowing
# down the automated testing

# TODO this is unused but maybe do not delete it
def split_to_full_splits(split, N):
    """
    Get all informative full splits compatible with the given partial split.
    @param split: a partial split of N leaves
    @param N: the number of leaves in the full tree
    @return: the set of all informative full splits compatible with the given split
    """
    accum = set()
    all_leaves = set(range(N))
    left, right = split
    free_leaves = all_leaves - flattened_split(split)
    # If the input split is a full split
    # then return the set that includes only that split.
    if not free_leaves:
        if is_informative(split):
            accum.add(split)
        return accum
    # Each free leaf can go to the left or right side.
    ordered_free_leaves = list(sorted(free_leaves))
    for i in range(2 ** len(free_leaves)):
        next_left = set(left)
        next_right = set(right)
        for j, leaf in enumerate(ordered_free_leaves):
            if (1<<j) & i:
                next_left.add(leaf)
            else:
                next_right.add(leaf)
        full_split = make_split((next_left, next_right))
        if is_informative(full_split):
            accum.add(full_split)
    return accum

def get_compatible_trees(split, available_trees):
    """
    @param split: a split is a frozenset of two frozensets of integers
    @param available_trees: a collection of trees
    @return: a subset of the available trees
    """
    compatible_subset = set()
    for tree in available_trees:
        if splitalg.split_tree_compatibility(split, tree):
            compatible_subset.add(tree)
    return compatible_subset

def get_splits_implied_by_tree(tree, available_splits):
    """
    This can be used to get the quartets implied by a tree.
    If an available partial split is implied by any
    full split of a tree, then the partial split will be in the returned set.
    """
    accum = set()
    for split in available_splits:
        if any(splitalg.split_implication(full_split, split) for full_split in tree):
            accum.add(split)
    return accum

def get_splits_implied_by_trees(trees, available_splits):
    """
    This can be used to get the quartets implied by a set of trees.
    If an available partial split is implied by any full split of each tree
    then the partial split will be in the returned set.
    """
    accum = set(available_splits)
    for tree in trees:
        accum = get_splits_implied_by_tree(tree, accum)
    return accum


class MultiHelper:
    """
    Use this for enumerating the multiconnected trees.
    Use it instead of passing massive numbers of arguments.
    """
    def __init__(self, N):
        self.N = N
        self.all_quartets = list(splitalg.get_quartets(N))
        self.all_splits = list(splitalg.get_informative_partial_splits(N))

def get_mc_trees_helper(current_split_set, compatible_trees, start_index, mh, mc_trees):
    """
    @param current_split_set: a set of less than N-3 splits
    @param compatible_trees: a set of trees that are compatible with the current set of splits
    @param start_index: search only (k choose 2) rather than (k*k).
    @param mh: a MultiHelper with data that is constant over the enumeration
    @param mc_trees: add to this list
    """
    nsplits = mh.N - 3
    for i in range(len(mh.all_splits) - start_index):
        target_split = mh.all_splits[start_index + i]
        if not any(splitalg.split_redundancy(target_split, split) for split in current_split_set):
            next_trees = get_compatible_trees(target_split, compatible_trees)
            if next_trees:
                next_set = current_split_set | set([target_split])
                # see if the set of compatible trees implies a quartet not implied by any member of the split set
                implied_quartets = get_splits_implied_by_trees(next_trees, mh.all_quartets)
                if splitalg.all_quartets_are_individually_implied(implied_quartets, next_set):
                    mc_trees.append(frozenset(next_set))
                if len(next_set) < nsplits:
                    get_mc_trees_helper(next_set, next_trees, start_index + i+1, mh, mc_trees)

def get_mc_trees_slow(N):
    """
    This is written much later in a way that I hope is clearer and slower.
    """
    if N <= 3:
        raise ValueError('not enough states')
    all_quartets = list(get_quartets(N))
    all_splits = list(get_all_partial_splits(N))
    mc_trees = []
    for nsplits in range(N-2):
        for splits in itertools.combinations(all_splits, nsplits):
            pass


def get_mc_trees(N):
    """
    A multiconnected tree is a frozenset of <= N-3 splits.
    The splits must satisfy various properties.
    The splits must be compatible.
    No split in the set implies another split in the set.
    The set of trees compatible with the split set must not
    imply a split not individually implied by a split in the set.
    @param N: the number of leaves
    """
    assert N > 3
    nsplits = N - 3
    # initialize the set of quartets and the set of informative splits given N leaves
    mh = MultiHelper(N)
    # get all possible bifurcating trees with N leaves
    trees = splitalg.get_bifurcating_trees(N)
    # try each combination of <= N-3 informative partial splits
    mc_trees = []
    initial_split_set = set()
    initial_compatible_trees = set(trees)
    initial_offset = 0
    get_mc_trees_helper(initial_split_set, initial_compatible_trees, initial_offset, mh, mc_trees)
    return mc_trees

def get_mc_trees_ben_helper(current_split_set, start_index, N, available_splits, mc_trees):
    """
    @param current_split_set: a set of less than N-3 splits
    @param start_index: search only (k choose 2) rather than (k*k).
    @param N: the number of leaves
    @param available_splits: a list of all possible informative splits of N leaves
    @param mc_trees: add to this list
    """

    # debug
    assert splitalg.pairwise_split_consistency(current_split_set)

    nsplits = N - 3
    for i in range(len(available_splits) - start_index):

        # debug
        assert splitalg.pairwise_split_consistency(current_split_set)

        target_split = available_splits[start_index + i]
        next_set = current_split_set | set([target_split])
        if all(splitalg.split_consistency(s, target_split) for s in current_split_set):

            # debug
            if not splitalg.pairwise_split_consistency(next_set):
                print 'FAIL'
                print 'current split set consistency:'
                print is_pairwise_rule_consistent(current_split_set)
                print 'current split set:'
                print current_split_set
                print 'target split:'
                print target_split
                print 'testing commutativity:'
                print all(splitalg.split_consistency(s, target_split) for s in current_split_set)
                print all(splitalg.split_consistency(target_split, s) for s in current_split_set)
                assert False

            mc_trees.append(frozenset(next_set))

            if len(next_set) < nsplits:

                # debug
                assert splitalg.pairwise_split_consistency(next_set)

                get_mc_trees_ben_helper(next_set, start_index + i+1, N, available_splits, mc_trees)

def get_mc_trees_ben(N):
    """
    Use Ben's nine rules to generate the set of mc trees with a given number of leaves.
    """
    assert N > 3
    nsplits = N - 3
    available_splits = list(splitalg.get_informative_partial_splits(N))
    mc_trees = []
    initial_split_set = set()
    initial_offset = 0
    get_mc_trees_ben_helper(initial_split_set, initial_offset, N, available_splits, mc_trees)
    return mc_trees

def get_mc_trees_ben_slow(N):
    """
    This is written much later in a way that I hope is clearer and slower.
    """
    if N <= 3:
        raise ValueError('not enough states')
    all_splits = list(splitalg.get_informative_partial_splits(N))
    mc_trees = []
    for nsplits in range(N-2):
        for splits in itertools.combinations(all_splits, nsplits):
            pass

def mc_tree_to_trees(mc_tree, available_trees):
    """
    Get the set of bifurcating trees implied by a multiconnected tree.
    @param mc_tree: a set of splits
    @param available_trees: all available bifurcating trees
    @return: a set of bifurcating trees
    """
    trees = available_trees
    for split in mc_tree:
        trees = get_compatible_trees(split, trees)
    return trees


class TestMCTree(unittest.TestCase):

    def foo_test_rules_trees(self):
        """
        Count mc trees computed using the nine rules.
        """
        print 'testing mc trees:'
        k = 5
        for i in range(4, 4+k):
            print 'N:', i
            all_bifurcating_trees = splitalg.get_bifurcating_trees(i)
            print len(all_bifurcating_trees), 'bifurcating trees'
            mc_trees_ben = set(get_mc_trees_ben(i))
            print len(mc_trees_ben), 'multiconnected trees'
            tree_set_sets_ben = set(frozenset(mc_tree_to_trees(tree, all_bifurcating_trees)) for tree in mc_trees_ben)
            print len(tree_set_sets_ben), 'induced sets of bifurcating trees'
        print

    def foo_test_counterexample(self):
        """
        Show an mc tree using the naive formulation but not the rules.
        """
        print 'counterexample commentary:'
        # see how the rules evaluates this set
        mc_tree = frozenset([
            frozenset([frozenset([1, 4]), frozenset([0, 2, 5])]),
            frozenset([frozenset([1, 2, 4, 5]), frozenset([0, 3])]),
            frozenset([frozenset([2, 3, 5]), frozenset([1, 4])])])
        a, b, c = mc_tree
        for alpha, beta in ((a, b), (b, c), (c, a)):
            result, comment = is_rule_consistent(alpha, beta, True)
            print comment
        # get the trees compatible with this set of splits
        N = 6
        trees = get_trees(N)
        for split in mc_tree:
            trees = get_compatible_trees(split, trees)
        print 'found', len(trees), 'compatible trees:'
        for tree in trees:
            print tree
        # find the quartets implied by the trees
        all_quartets = get_quartets(N)
        quartets = get_splits_implied_by_trees(trees, all_quartets)
        print len(quartets), 'quartets implied by the set of bifurcating trees:'
        for q in quartets:
            print q
        print

    def foo_test_counterexample_b(self):
        print 'test this counterexample'
        A = make_split(([0,1], [2,3]))
        B = make_split(([0,1], [4,5]))
        invalid_tree = frozenset([A, B])
        # get the trees compatible with this set of splits
        N = 6
        trees = get_trees(N)
        for split in invalid_tree:
            trees = get_compatible_trees(split, trees)
        print len(trees), 'compatible trees:'
        for tree in trees:
            print tree
        # find the quartets implied by the trees
        all_quartets = get_quartets(N)
        quartets = get_splits_implied_by_trees(trees, all_quartets)
        print len(quartets), 'quartets implied by the set of bifurcating trees:'
        for q in quartets:
            print q
        print

    def foo_test_mc_trees(self):
        print 'testing mc trees:'
        k = 4
        for i in range(4, 4+k):
            print 'N:', i
            all_bifurcating_trees = splitalg.get_bifurcating_trees(i)
            all_quartets = splitalg.get_quartets(i)
            # get multifurcating trees using the nine rules
            mc_trees_ben = set(get_mc_trees_ben(i))
            print 'number of mc trees using the nine rules:', len(mc_trees_ben)
            tree_set_sets_ben = set(frozenset(mc_tree_to_trees(tree, all_bifurcating_trees)) for tree in mc_trees_ben)
            print 'number of unique bifurcating tree sets generated by these mc trees:', len(tree_set_sets_ben)
            # get multifurcating trees using a naive method
            mc_trees = set(get_mc_trees(i))
            print 'number of mc trees using a naive search:', len(mc_trees)
            tree_set_sets = set(frozenset(mc_tree_to_trees(tree, all_bifurcating_trees)) for tree in mc_trees)
            print 'number of unique bifurcating tree sets generated by these mc trees:', len(tree_set_sets)
            # compare subsets
            symmetric_difference = mc_trees_ben ^ mc_trees
            print 'number of trees in the symmetric difference set:', len(symmetric_difference)
            if mc_trees_ben - mc_trees:
                print 'an mc tree using the rules but not the naive:'
                print (mc_trees_ben - mc_trees).pop()
            if mc_trees - mc_trees_ben:
                print 'an mc tree using the naive but not the rules:'
                print (mc_trees - mc_trees_ben).pop()
            # Look for a naive mc tree that generates a set of bifurcating trees that
            # is not representable by a multiconnected tree constructed from the rules.
            weird_naive_mc_trees = set()
            for mc_tree in mc_trees:
                if frozenset(mc_tree_to_trees(mc_tree, all_bifurcating_trees)) not in tree_set_sets_ben:
                    weird_naive_mc_trees.add(mc_tree)
            if weird_naive_mc_trees:
                weird_tree = weird_naive_mc_trees.pop()
                print 'a naive mc tree that represents a set of trees not representable by any real mc tree:'
                print weird_tree
                mytrees = set(all_bifurcating_trees)
                for split in weird_tree:
                    mytrees = get_compatible_trees(split, mytrees)
                print len(mytrees), 'trees compatible with the weird tree:'
                for tree in mytrees:
                    print tree
                # find the quartets implied by the trees
                quartets = get_splits_implied_by_trees(mytrees, all_quartets)
                print len(quartets), 'quartets implied by the set of bifurcating trees:'
                for q in quartets:
                    print q
            print

    def foo_test_counterexample_e(self):
        """
        Show an mc tree using the naive formulation but not the rules.
        """
        print 'counterexample commentary:'
        # see how the rules evaluates this set
        mc_tree = frozenset([
            frozenset([frozenset([1, 4]), frozenset([0, 2, 5])]),
            frozenset([frozenset([1, 2, 4, 5]), frozenset([0, 3])]),
            frozenset([frozenset([2, 3, 5]), frozenset([1, 4])])])
        a, b, c = mc_tree
        for alpha, beta in ((a, b), (b, c), (c, a)):
            result, comment = split_consistency(alpha, beta, verbose=True)
            print comment
        # get the trees compatible with this set of splits
        N = 6
        trees = get_bifurcating_trees(N)
        for split in mc_tree:
            trees = set(t for t in trees if split_tree_compatibility(split, t))
        print 'found', len(trees), 'compatible trees:'
        for tree in trees:
            print tree
        # find the quartets implied by the trees
        qs = get_quartets(N)
        for t in trees:
            quartets = set(q for q in qs if split_tree_compatibility(q, t))
        print len(qs), 'quartets implied by the set of bifurcating trees:'
        for q in qs:
            print q
        print 'splits defining the mc-tree:'
        print mc_tree
        print 'all quartets individually implied:',
        print all_quartets_are_individually_implied(qs, mc_tree)


if __name__ == '__main__':
    unittest.main()
