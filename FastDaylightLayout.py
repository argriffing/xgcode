"""
Use a C extension module to do fast tree layout.
"""

import unittest
from optparse import OptionParser

import SpatialTree
import EqualArcLayout
import day

class LayoutError(Exception): pass


def _build_dtree(dtree, node, count):
    """
    @param dtree: an object that lives in my C extension
    @param node: like a SpatialNode
    @param count: the number of nodes added so far
    """
    x, y = node.location
    node.dtree_id = count
    count += 1
    dtree.begin_node(node.dtree_id)
    dtree.set_x(x)
    dtree.set_y(y)
    for child in node.children:
        count = _build_dtree(dtree, child, count)
    dtree.end_node()
    return count

def get_neighbor_count(node):
    neighbor_count = len(node.children)
    if node.parent:
        neighbor_count += 1
    return neighbor_count


class StraightBranchLayout:

    def __init__(self):
        self.iteration_count = 3

    def set_iteration_count(self, iteration_count):
        """
        @param iteration_count: do this many daylight equalizing iterations
        """
        self.iteration_count = iteration_count

    def do_layout(self, tree):
        """
        @param tree: something like a SpatialTree
        """
        # create the initial layout
        EqualArcLayout.do_layout(tree)
        # use sax-like events to create a parallel tree in the C extension
        dtree = day.Day()
        count = _build_dtree(dtree, tree.root, 0)
        # repeatedly reroot and equalize
        for iteration in range(self.iteration_count):
            for node in tree.breadth_first():
                neighbor_count = len(node.children)
                if node.parent:
                    neighbor_count += 1
                if neighbor_count > 2:
                    dtree.select_node(node.dtree_id)
                    dtree.reroot()
                    try:
                        dtree.equalize()
                    except RuntimeError as e:
                        raise LayoutError(e)
        # extract the x and y coordinates from the parallel tree
        for node in tree.preorder():
            dtree.select_node(node.dtree_id)
            x = dtree.get_x()
            y = dtree.get_y()
            node.location = (x, y)
        # take off the silly dtree_id members
        for node in tree.preorder():
            del node.dtree_id


class CurvedBranchLayout:

    def __init__(self):
        self.min_segment_count = 200
        self.node_factory = SpatialTree.SpatialTreeNode

    def set_min_segment_count(self, min_segment_count):
        """
        The tree will be automatically segmented.
        @param min_segment_count: the minimum number of segments
        """
        self.min_segment_count = min_segment_count

    def set_node_factory(self, node_factory):
        """
        @param node_factory: a way to generate nodes that segment the branches
        """
        self.node_factory = node_factory

    def do_layout(self, tree):
        """
        Interleave the processes of breaking the tree and daylight equalization.
        @param tree: something like a SpatialTree
        """
        # reroot to a node with more than two neighbors
        if get_neighbor_count(tree.root) < 3:
            suitable_roots = [
                    x for x in tree.preorder() if get_neighbor_count(x) > 2]
            if suitable_roots:
                tree.reroot(suitable_roots[0])
        # create the initial layout
        EqualArcLayout.do_layout(tree)
        # determine the maximum allowed branch length
        total_branch_length = tree.get_total_length()
        max_branch_length = total_branch_length / float(self.min_segment_count)
        # any branch longer than the max branch length will be broken in half
        iteration = 0
        while True:
            print 'progressive iteration', iteration+1
            # write the extension tree
            dtree = day.Day()
            _build_dtree(dtree, tree.root, 0)
            # equalize the layout
            for node in tree.breadth_first():
                if get_neighbor_count(node) > 1:
                    dtree.select_node(node.dtree_id)
                    dtree.reroot()
                    try:
                        dtree.equalize()
                    except RuntimeError as e:
                        raise LayoutError(e)
            # read the extension tree
            for node in tree.preorder():
                dtree.select_node(node.dtree_id)
                x = dtree.get_x()
                y = dtree.get_y()
                node.location = (x, y)
            # break the branches
            old_nodes = list(tree.preorder())
            for node in old_nodes:
                if node is tree.root:
                    if node.blen is not None:
                        raise HandlingError(
                                'the root node should not have a branch length')
                elif node.blen is None:
                    raise HandlingError(
                            'each non-root node should have a branch length')
                elif node.blen > max_branch_length:
                    # create a new node and set its attributes
                    new = self.node_factory()
                    new.name = node.name
                    # insert the new node
                    tree.insert_node(new, node.parent, node, .5)
            # if no node was added then break out of the loop
            if len(old_nodes) == len(list(tree.preorder())):
                break
            else:
                iteration += 1



def segment_tree(tree, min_segment_count, node_factory):
    """
    Break a tree into some minimum number of segments.
    @param tree: a tree with branch lengths
    @param min_segment_count: the minimum number of segments for the tree
    """
    # determine the maximum allowed branch length
    total_branch_length = tree.get_total_length()
    max_branch_length = total_branch_length / float(min_segment_count)
    # any branch longer than the max branch length will be broken in half
    while True:
        old_nodes = list(tree.preorder())
        for node in old_nodes:
            if node is tree.root:
                if node.blen is not None:
                    raise HandlingError(
                            'the root node should not have a branch length')
            elif node.blen is None:
                raise HandlingError(
                        'each non-root node should have a branch length')
            elif node.blen > max_branch_length:
                # create a new node and set its attributes
                new = node_factory()
                new.name = node.name
                # insert the new node
                tree.insert_node(new, node.parent, node, .5)
        # if no node was added then break out of the loop
        if len(old_nodes) == len(list(tree.preorder())):
            break

def make_file(newick_string, filename_prefix, iterations, min_segments):
    import DrawTreeImage
    import Newick
    import SpatialTree
    # make a spatial tree
    tree = Newick.parse(newick_string, SpatialTree.SpatialTree)
    # break the tree into pieces
    segment_tree(tree, min_segments, SpatialTree.SpatialTreeNode)
    print 'broke the tree into', sum(1 for node in tree.preorder()), 'nodes'
    # do the layout
    do_layout(tree, iterations)
    print 'did the layout'
    # color some of the branches
    red = 'ff0000'
    green = '00ff00'
    blue = '0000ff'
    for child, color in zip(tree.root.children, (red, green, blue)):
        for node in child.preorder():
            node.branch_color = color
    # get the image string
    image_format = 'png'
    image_string = DrawTreeImage.get_tree_image(tree, (640, 480), image_format)
    print 'created the image string'
    # write the image file
    with open('%s.%s' % (filename_prefix, image_format), 'wb') as fout:
        fout.write(image_string)

def make_file_progressive(newick_string, filename_prefix, min_segment_count):
    import DrawTreeImage
    import Newick
    import SpatialTree
    # make a spatial tree
    tree = Newick.parse(newick_string, SpatialTree.SpatialTree)
    # do the layout
    layout = CurvedBranchLayout()
    layout.set_min_segment_count(min_segment_count)
    layout.do_layout(tree)
    # color some of the branches
    red = 'ff0000'
    green = '00ff00'
    blue = '0000ff'
    for child, color in zip(tree.root.children, (red, green, blue)):
        for node in child.preorder():
            node.branch_color = color
    # get the image string
    image_format = 'png'
    image_string = DrawTreeImage.get_tree_image(tree, (640, 480), image_format)
    print 'created the image string'
    # write the image file
    with open('%s.%s' % (filename_prefix, image_format), 'wb') as fout:
        fout.write(image_string)

def make_hiv_images():
    import urllib
    url_name_pairs = [
            ('http://www.bioinf.manchester.ac.uk/ctree/sampleData/gagRefTree.ph', 'gag'),
            ('http://www.bioinf.manchester.ac.uk/ctree/sampleData/polRefTree.ph', 'pol'),
            ('http://www.bioinf.manchester.ac.uk/ctree/sampleData/envRefTree.ph', 'env')
            ]
    for url, name in url_name_pairs:
        fin = urllib.urlopen(url)
        newick_string = fin.read()
        fin.close()
        make_file_progressive(newick_string, name, 500)


class TestLayout(unittest.TestCase):

    def test_progressive(self):
        branch_manager_string = "(((((((((PV22:0.3,BH10:0.3):0.1,BRU:0.5):0.1,HXB:0.7):2.4,SF2:3.3):0.1,CDC:3.7):0.5,WMJ2:3.0):0.4,RF:4.3):2.6,(ELI:6.3,MAL:6.1):1.9):2.7,Z3:9.3);"
        branch_manager_string_tri = "((((((((PV22:0.3,BH10:0.3):0.1,BRU:0.5):0.1,HXB:0.7):2.4,SF2:3.3):0.1,CDC:3.7):0.5,WMJ2:3.0):0.4,RF:4.3):2.6,(ELI:6.3,MAL:6.1):4.6,Z3:9.3);"
        make_file_progressive(branch_manager_string, 'branch-manager', 500)


if __name__ == '__main__':
    unittest.main()
