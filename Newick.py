"""
Represent newick trees.
"""

from optparse import OptionParser
import unittest

import Util

# base class for NewickTree
import Tree

# example trees
from NewickIO import rooted_example_tree
from NewickIO import daylight_example_tree
from NewickIO import brown_example_tree
from NewickIO import ben_example_tree

# tree writers
from NewickIO import get_newick_string
from NewickIO import get_multiline_newick_string
from NewickIO import get_narrow_newick_string

# tree readers
from NewickIO import parse

# tree reader exceptions
from NewickIO import NewickSyntaxError

# a generic test case
from NewickIO import TestNewick


class NewickSearchError(Exception):
    pass


class NewickNode(Tree.TreeNode):
    """
    Each node in a newick string has a name, zero to three children, and possibly a branch length.
    """

    def __init__(self):
        Tree.TreeNode.__init__(self)
        self.name = None
        self.blen = None

    def get_name(self):
        return self.name

    def set_branch_length(self, branch_length):
        self.blen = branch_length

    def get_branch_length(self):
        return self.blen

    def set_name(self, name):
        #TODO cf add_name
        if self.name is not None:
            raise NewickSyntaxError('the node already has a name')
        self.name = name

    def add_name(self, name):
        #TODO cf set_name
        if self.name is not None:
            raise NewickSyntaxError('each newick node must have only one name')
        self.name = name

    def __str__(self):
        if self.get_name() is not None:
            return str(self.get_name())
        else:
            return 'node_%d' % id(self)


class NewickTree(Tree.Tree):
    """
    This class represents a naive newick tree.
    """

    NodeFactory = NewickNode

    def __init__(self, root=None):
        Tree.Tree.__init__(self, root)

    def get_newick_string(self):
        return get_newick_string(self)

    def get_total_length(self):
        """
        Return the sum of all branch lengths on the tree.
        """
        return sum(node.blen for node in self.preorder() if node.blen)

    def gen_description_lines(self):
        """
        Yield single line strings describing features of the tree.
        """
        for description_line in Tree.Tree.gen_description_lines(self):
            yield description_line
        if self.has_branch_lengths():
            nbranches_with_lengths = len([node for node in self.gen_non_root_nodes() if node.blen is not None])
            yield 'the number of branches for which a length is available is %d' % nbranches_with_lengths
        else:
            nbranches_without_lengths = len([node for node in self.gen_non_root_nodes() if node.blen is None])
            yield 'the number of branches for which a length is not available is %d' % nbranches_without_lengths
        yield 'the number of branches with zero length is %d' % len([node for node in self.gen_non_root_nodes() if node.blen == 0.0])
        yield 'the diameter of the tree is %f' % self.get_diameter()


    def get_unique_node(self, name):
        valid_nodes = [node for node in self.preorder() if node.name == name]
        if not valid_nodes:
            raise NewickSearchError('no node named "%s" was found' % name)
        elif len(valid_nodes) > 1:
            raise NewickSearchError('%d nodes named "%s" were found' % (len(valid_nodes), name))
        return valid_nodes[0]

    def assert_valid(self):
        """
        Assert that the tree is a valid newick tree.
        These extra assertions are not required by the standard newick definition.
        """
        if not self.root:
            raise NewickSyntaxError('the newick tree must have at least one node')
        if self.has_branch_lengths():
            for child in self.root.children:
                for node in child.preorder():
                    if node.blen is None:
                        raise NewickSyntaxError('if any node has a branch length then all non-root nodes must have a branch length')

    def has_negative_branch_lengths(self):
        for node in self.gen_non_root_nodes():
            if node.blen is not None:
                if node.blen < 0:
                    return True
        return False

    def add_branch_lengths(self):
        """
        Add branch lengths if the tree does not have any.
        """
        if not self.has_branch_lengths():
            for node in self.gen_non_root_nodes():
                node.blen = 1.0

    def has_branch_lengths(self):
        for node in self.preorder():
            if node.blen is not None:
                return True
        return False

    def is_rooted(self):
        if not self.root:
            return False
        if len(self.root.children) < 3:
            return True
        return False

    def prune(self, node):
        """
        Remove a tip of the tree and prune the tree all the way back to an intersection.
        If the root is on the path back to an intersection then the root will be moved.
        @param node: the node that defines the tip of the branch to be pruned
        @return: the node that is the intersection with the rest of the tree
        """
        while node.get_neighbor_count() == 1:
            if node.parent:
                parent = node.parent
                parent.children = [child for child in parent.children if child is not node]
                node.parent = None
                node = parent
            else:
                child = node.children[0]
                node.children = []
                child.parent = None
                child.blen = None
                self.root = child
                node = child
        return node

    def remove_node(self, node):
        """
        Add the branch lengths of the node to its child nodes and move the child nodes up a level.
        Removing the root node is not allowed.
        The complexity of this operation is independent of the size of the tree.
        @param node: the node to be removed
        """
        parent = node.parent
        if not parent:
            if len(node.children) == 1:
                self.root = node.children[0]
            return
        # add the branch length of the node to each child node
        if node.blen is not None:
            for child in node.children:
                child.blen += node.blen
        # clear the branch length of the removed node
        node.blen = None
        # remove the node from its parent
        parent.children = [child for child in parent.children if child is not node]
        node.parent = None
        # attach the children to the parent
        for child in node.children:
            parent.children.append(child)
            child.parent = parent

    def insert_node(self, node, parent, child, fraction):
        """
        Insert a new node between the parent and child nodes.
        This operation is useful for rerooting the tree.
        The complexity of this operation is independent of the size of the tree.
        @param node: the new node to be inserted
        @param parent: the parent of the new node
        @param child: the child of the new node
        @param fraction: a number between 0 and 1 that is small when the inserted node is closer to the parent
        """
        assert child in parent.children
        assert 0 <= fraction <= 1
        # update the branch lengths
        if child.blen is not None:
            node.blen = fraction * child.blen
            child.blen = (1 - fraction) * child.blen
        # update the parent node
        siblings = parent.children[:]
        parent.children = []
        for sibling in siblings:
            if sibling is child:
                parent.children.append(node)
            else:
                parent.children.append(sibling)
        # update the inserted node
        node.children = [child]
        node.parent = parent
        # update the child node
        child.parent = node

    def reroot(self, node):
        """
        Reorder the tree so that the specified preexisting node is at the root.
        This function does not change the number of nodes in the tree.
        @param node: this node should already exist within the tree
        """
        assert node
        if node is self.root:
            return
        nodes = list(reversed(self.get_path_to_root(node)))
        for parent, child in zip(nodes[:-1], nodes[1:]):
            child.blen, parent.blen = parent.blen, child.blen
            parent.children.remove(child)
            child.children.append(parent)
            parent.parent = child
            child.parent = None
        self.root = node

    def get_diameter(self):
        """
        @return: the maximum pairwise distance between leaves
        """
        # TODO add a test for this function
        # get the heights
        id_to_height = {}
        for node in self.postorder():
            # initialize height to zero
            height = 0
            # if the node has children then add the height of the tallest child
            if node.has_children():
                height = max(id_to_height[id(child)] for child in node.gen_children())
            # if the node has a branch length then add the branch length
            blen = node.get_branch_length()
            if blen:
                height += blen
            # cache the height
            id_to_height[id(node)] = height
        # get the diameters
        id_to_diameter = {}
        for node in self.postorder():
            candidates = [0]
            # one candidate for the diameter is the sum of the heights of the two tallest subtrees
            if node.get_child_count() > 1:
                subtree_heights = [id_to_height[id(child)] for child in node.gen_children()]
                negative_heights = [-height for height in subtree_heights]
                candidate_height = -(Util.select(negative_heights, 0) +  Util.select(negative_heights, 1))
                candidates.append(candidate_height)
            # other candidates are the diameters of the subtrees
            for child in node.gen_children():
                candidates.append(id_to_diameter[id(child)])
            # set the diameter
            id_to_diameter[id(node)] = max(candidates)
        # return the diameter associated with the root
        return id_to_diameter[id(self.get_root())]


class TestNewickDerived(TestNewick):
    TreeFactory = NewickTree


def main():
    import DrawTree
    print 'input newick string:'
    unrooted_example_tree = daylight_example_tree
    print unrooted_example_tree
    tree = parse(unrooted_example_tree, NewickTree)
    tree.assert_valid()
    print 'ascii tree representation:'
    drawer = DrawTree.DrawTree()
    drawer.use_branch_lengths = True
    drawer.vertical_spacing = 1
    drawer.horizontal_spacing = 2
    print drawer.draw(tree)
    print '\n'.join(tree.gen_description_lines())
    print 'output newick string:'
    print tree.get_newick_string()


if __name__ == '__main__':
    parser = OptionParser()
    #parser.add_option('-v', '--verbose', action='store_true', dest='verbose', default=False)
    #parser.add_option('-o', '--output', dest='output_filename', metavar='FILE', help='output file')
    parser.add_option('--test', action='store_true', dest='test', default=False)
    options, args = parser.parse_args()
    if options.test:
        suite = unittest.TestLoader().loadTestsFromTestCase(TestNewickDerived)
        unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        main()

