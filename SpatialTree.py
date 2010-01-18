#!/usr/bin/env python

import math
import unittest
from optparse import OptionParser

import Newick

# a utility variable defining 360 degrees
m2pi = 2*math.pi


def get_angle(pa, pb):
    """
    @param pa: an (x, y) pair
    @param pb: an (x, y) pair
    @return: the angle towards point pb from point pa
    """
    ax, ay = pa
    bx, by = pb
    return math.atan2(by - ay, bx - ax)


class AngleIntervalError(Exception):
    pass


class AngleInterval:
    """
    An angle interval is a (low, high) pair.
    The low and high values are each in the real number interval [0, 2*pi).
    """

    def __init__(self, low, high):
        self.low = low % m2pi
        self.high = high % m2pi

    def get_magnitude(self):
        """
        @return: the angle spanned by the interval
        """
        return (self.high - self.low) % m2pi

    def update(self, other):
        """
        Modify the current interval by adding another interval.
        The union of the ranges is assumed to be contiguous and to span less than 2*pi radians.
        """
        triples = []
        for low in (self.low, other.low):
            for high in (self.high, other.high):
                triples.append(((high-low) % m2pi, low, high))
        best_magnitude, best_low, best_high = max(triples)
        self.low, self.high = best_low, best_high

    def contains_angle(self, theta):
        if self.low == self.high:
            return False
        phi = 0
        phi += (self.high - theta) % m2pi
        phi += (theta - self.low) % m2pi
        return phi < m2pi



def _get_scaling_factor(current_size, max_size):
    """
    @param current_size: the width and height of the bounding box of the image we want to draw
    @param max_size: the maximum allowed width and height
    @return: the amount by which the current size should be scaled
    """
    cwidth, cheight = current_size
    mwidth, mheight = max_size
    xscale = float(mwidth) / float(cwidth)
    yscale = float(mheight) / float(cheight)
    return min(xscale, yscale)


class SpatialTreeBranch:
    def __init__(self, src_location, dst_location, node):
        """
        @param src_location: an (x, y) location on the display
        @param dst_location: an (x, y) location on the display
        @param node: the node associated with the branch
        """
        self.src_location = src_location
        self.dst_location = dst_location
        self.branch_color = getattr(node, 'branch_color', None)


class SpatialTreeNode(Newick.NewickNode):

    def __init__(self):
        Newick.NewickNode.__init__(self)
        # theta is relative to the direction of the parent node
        # this is obsolescent
        self.theta = None
        # location is with respect to the layout coordinate system not the display coordinate system
        self.location = None

    def set_location(self, location):
        """
        @param location: an (x, y) coordinate pair that defines a point
        """
        self.location = location

    def get_location(self):
        """
        @return: an (x, y) coordinate pair that defines a point
        """
        return self.location


class SpatialTree(Newick.NewickTree):
    """
    A spatial tree has information about the angles between branches.
    """

    NodeFactory = SpatialTreeNode

    def __init__(self, root=None):
        """
        The 2d transformation referred to here is for putting the tree into a box.
        """
        Newick.NewickTree.__init__(self, root)
        self.scale = 1
        self.theta = 0
        self.center = (0, 0)

    def _layout_to_display(self, location):
        ct = math.cos(self.theta)
        st = math.sin(self.theta)
        cx, cy = self.center
        x, y = location
        x, y = x*ct - y*st, y*ct + x*st
        x, y = self.scale*x, self.scale*y
        x, y = x - cx, y-cy
        return (x, y)

    def fit(self, max_size):
        """
        Set the display transformation parameters to fit in a box of max_size.
        Search angles in a range of 180 degrees.
        This function sets a new scale, a new rotation, and a new translation.
        """
        self.scale = 1
        self.theta = 0
        self.center = (0, 0)
        sf_theta_pairs = []
        increment_count = 60
        for i in range(increment_count):
            self.theta = i * math.pi / increment_count
            xmin, ymin, xmax, ymax  = self.get_extents()
            cx = xmax - xmin
            cy = ymax - ymin
            if cx == 0:
                cx = cy / 100
            if cy == 0:
                cy = cx / 100
            current_size = (cx, cy)
            sf = _get_scaling_factor(current_size, max_size)
            sf_theta_pairs.append((sf, self.theta))
        self.scale, self.theta = max(sf_theta_pairs)
        xmin, ymin, xmax, ymax  = self.get_extents()
        self.center = ((xmin+xmax)/2.0, (ymin+ymax)/2.0)

    def get_extents(self):
        """
        This call assumes that a layout has already been done.
        @return: (min(x), min(y), max(x), max(y))
        """
        tip_locations = [self._layout_to_display(tip.location) for tip in self.gen_tips()]
        xs = [x for x, y in tip_locations]
        ys = [y for x, y in tip_locations]
        return (min(xs), min(ys), max(xs), max(ys))

    def gen_branches(self):
        """
        This call assumes that a layout has already been done.
        Each node has a corresponding branch except for the root node.
        """
        for node in self.gen_non_root_nodes():
            src_display_location = self._layout_to_display(node.parent.location)
            dst_display_location = self._layout_to_display(node.location)
            yield SpatialTreeBranch(src_display_location, dst_display_location, node)

    def insert_node(self, node, parent, child, fraction):
        """
        Insert a new node between the parent and child nodes.
        @param node: the new node to be inserted
        @param parent: the parent of the new node
        @param child: the child of the new node
        @param fraction: a number between 0 and 1 that is small when the inserted node is closer to the parent
        """
        Newick.NewickTree.insert_node(self, node, parent, child, fraction)
        if parent.location and child.location:
            parent_x, parent_y = parent.location
            child_x, child_y = child.location
            node_x = fraction * child_x + (1 - fraction) * parent_x
            node_y = fraction * child_y + (1 - fraction) * parent_y
            node.location = (node_x, node_y)


class TestSpatialTree(unittest.TestCase):

    def test_basic_parsing(self):
        tree = Newick.parse(Newick.daylight_example_tree, SpatialTree)

    def test_degenerate_angle_interval(self):
        angle_interval = AngleInterval(0, 0)
        self.assertFalse(angle_interval.contains_angle(0))
        self.assertFalse(angle_interval.contains_angle(1.234))

    def test_small_angle_interval(self):
        angle_interval = AngleInterval(5.5, 1.1)
        self.assertTrue(angle_interval.contains_angle(0))
        self.assertFalse(angle_interval.contains_angle(2.2))
        self.assertFalse(angle_interval.contains_angle(4.4))


def main():
    import EqualArcLayout
    import DrawTreeImage
    # read a tree
    tree_string = Newick.daylight_example_tree
    tree = Newick.parse(tree_string, SpatialTree)
    # do the layout
    EqualArcLayout.do_layout(tree)
    # show some debugging information
    max_size = (640, 480)
    tree.fit(max_size)
    extents = tree.get_extents()
    print 'extents:', extents
    branches = list(tree.gen_branches())
    print 'branches:'
    for branch in branches:
        print branch.src_location, branch.dst_location
    # draw the image
    image_format = 'png'
    image_string = DrawTreeImage.get_tree_image(tree, max_size, image_format)
    # write the image file
    fout = open('test.%s' % image_format, 'wb')
    fout.write(image_string)
    fout.close()

if __name__ == '__main__':
    parser = OptionParser()
    #parser.add_option('-v', '--verbose', action='store_true', dest='verbose', default=False)
    #parser.add_option('-o', '--output', dest='output_filename', metavar='FILE', help='output file')
    parser.add_option('--test', action='store_true', dest='test', default=False)
    options, args = parser.parse_args()
    if options.test:
        suite = unittest.TestLoader().loadTestsFromTestCase(TestSpatialTree)
        unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        main()



