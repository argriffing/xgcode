"""
Simulate some unrooted bifurcating trees.
"""

import unittest
import random
import math

import numpy as np

import NewickIO
import FelTree
import Xtree

def sample_agglomerated_tree(ntaxa):
    """
    Sample a weighted phylogenetic xtree.
    All of the branch lengths will be default lengths.
    This method follows the simulation method used in "Why neighbor-joining works".
    It agglomerates subtrees at random.
    @param ntaxa: the number of leaves in the tree
    @return: the root of a weighted phylogenetic xtree
    """
    # a tree must have at least three leaves before it has an internal node
    assert ntaxa > 2
    # initialize the pool of subtrees
    subtrees = []
    for i in range(ntaxa):
        v = Xtree.WPXVertex()
        v.label = i
        subtrees.append(v)
    # repeatedly agglomerate pairs of subtrees
    while len(subtrees) > 3:
        root = Xtree.WPXVertex()
        # select the items and efficiently delete them from the list
        for i in range(2):
            a = random.randrange(len(subtrees))
            root.add_child(subtrees[a])
            subtrees[a] = subtrees[-1]
            del subtrees[-1]
        # add the new subtree to the list
        subtrees.append(root)
    # agglomerate the final three subtrees
    root = Xtree.WPXVertex()
    for t in subtrees:
        root.add_child(t)
    return root

def sample_tree(n_base_leaves, n_expected_extra_leaves, expected_branch_length):
    # calculate parameters that give the requested expected values
    poisson_lambda = n_expected_extra_leaves
    gamma_shape = math.sqrt(expected_branch_length)
    gamma_scale = math.sqrt(expected_branch_length)
    # get the number of leaves by rejection sampling
    while True:
        nleaves = n_base_leaves + np.random.poisson(poisson_lambda)
        if (3 <= nleaves <= 52):
            break
    # create a bunch of unconnected nodes
    pool = []
    for i in range(nleaves):
        # make a leaf
        leaf = FelTree.NewickNode()
        # name the leaf
        name = chr(ord('a') + (i % 26))
        if i >= 26:
            name = name.upper()
        leaf.add_name(name)
        # add the leaf to the pool
        pool.append(leaf)
    # repeatedly take two subtrees from the pool and replace them with a single subtree
    while len(pool) > 1:
        root = FelTree.NewickNode()
        for i in range(2):
            child = pool.pop(random.randrange(len(pool)))
            root.add_child(child)
        pool.append(root)
    # the tree should be unrooted instead of rooted
    tree = FelTree.NewickTree(pool[0])
    new_root = max((node.get_neighbor_count(), node) for node in tree.preorder())[1]
    old_root = tree.get_root()
    tree.set_root(new_root)
    tree.remove_node(old_root)
    # sample branch lengths
    for node in tree.preorder():
        if not node.is_root():
            branch_length = np.random.gamma(gamma_shape, gamma_scale)
            node.set_branch_length(branch_length)
    # return the tree
    return tree


class TestTreeSampler(unittest.TestCase):

    def test_tree_sampler(self):
        n_base_leaves = 4
        n_expected_extra_leaves = 1
        expected_branch_length = 1
        tree = sample_tree(n_base_leaves, n_expected_extra_leaves, expected_branch_length)

    def test_xtree_sampler(self):
        ntaxa = 20
        root = sample_agglomerated_tree(ntaxa)
        self.assertEquals(len(root.get_labeled_vertices()), ntaxa)
        self.assertEquals(len(root.get_branches()), ntaxa*2 - 3)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestTreeSampler)
    unittest.TextTestRunner(verbosity=2).run(suite)
