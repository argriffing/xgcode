"""
Do things with the Kingman coalescent.

This module uses similar notation to the Ftree module.
The following blurb is from a random slide on the internet.
Kingman coalescent process for constructing a tree of k gene copies.
Go back in an amount of time t chosen from an exponential 4N/(k(k-1)).
The units of t is generations.
Combine two randomly chosen lineages.
Decrease k by one.
Iterate until there is just one left.
"""

import unittest
import random

import Ftree
import iterutils


def sample(k):
    """
    Branch lengths of the returned tree have strange units.
    They are in generations per population size.
    So multiplying the branch lengths by the population size
    will give branch lengths in units of generations.
    @param k: the number of leaves
    @return: (R, B) using Ftree notation
    """
    t = 0
    roots = set(range(k))
    r_next = k
    v_to_age = dict((r, t) for r in roots)
    R = set()
    B = {}
    while len(roots) > 1:
        # define the coalescence rate
        k = len(roots)
        coalescence_rate = k * (k - 1)
        # sample the next coalescence time
        t_delta = random.expovariate(coalescence_rate)
        t += t_delta
        # sample the pair of coalescing lineages
        va = random.choice(list(roots))
        roots.remove(va)
        vb = random.choice(list(roots))
        roots.remove(vb)
        # update the R and B structures
        r = r_next
        r_next += 1
        v_to_age[r] = t
        roots.add(r)
        R.add((r, va))
        R.add((r, vb))
        B[frozenset([r, va])] = t - v_to_age[va]
        B[frozenset([r, vb])] = t - v_to_age[vb]
    return R, B

def get_paths_to_root(R):
    sources, sinks = zip(*R)
    leaves = set(sinks) - set(sources)
    v_to_source = Ftree.R_to_v_to_source(R)
    paths_to_root = []
    for v in leaves:
        path = [v]
        while path[-1] in v_to_source:
            path.append(v_to_source[path[-1]])
        paths_to_root.append(path)
    return paths_to_root

def RB_to_v_to_age(R, B):
    """
    @param R: directed topology
    @param B: branch lengths in time units
    @return: map from vertex to age
    """
    sources, sinks = zip(*R)
    leaves = set(sinks) - set(sources)
    v_to_age = dict((v, 0) for v in leaves)
    v_to_source = Ftree.R_to_v_to_source(R)
    for v in Ftree.R_to_postorder(R):
        p = v_to_source.get(v, None)
        if p is not None:
            v_to_age[p] = v_to_age[v] + B[frozenset([v, p])]
    return v_to_age


class TestKingman(unittest.TestCase):

    def test_nbranches(self):
        nleaves = 10
        nbranches = 2 * nleaves - 2
        R, B = sample(nleaves)
        self.assertEqual(len(R), nbranches)
        self.assertEqual(len(B), nbranches)

    def test_clocklike(self):
        nleaves = 10
        R, B = sample(nleaves)
        paths_to_root = get_paths_to_root(R)
        ages = []
        for path in paths_to_root:
            age = sum(B[frozenset(p)] for p in iterutils.pairwise(path))
            ages.append(age)
        self.assertEqual(len(ages), nleaves)
        self.assertEqual(len(set(ages)), 1)


if __name__ == '__main__':
    unittest.main()

