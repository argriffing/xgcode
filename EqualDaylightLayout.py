"""
This is an implementation of the equal daylight layout algorithm.
The algorithm is described on pages 582-584 of Inferring Phylogenies.
This module should only be accessed through its do_layout() function.
"""

from optparse import OptionParser
import unittest
import math
import random
import profile

import cairo

import EqualArcLayout
import SpatialTree


# a utility variable defining 360 degrees
m2pi = 2*math.pi


def get_subtree_angle_interval(spatial_tree_node, center):
    # TODO figure out what these parameter values are exactly and add a docstring
    if center == spatial_tree_node.get_parent().get_location():
        angle = SpatialTree.get_angle(center, spatial_tree_node.get_location())
        interval = SpatialTree.AngleInterval(angle, angle)
    else:
        aa = SpatialTree.get_angle(center, spatial_tree_node.get_location())
        ab = SpatialTree.get_angle(center, spatial_tree_node.get_parent().get_location())
        if (ab - aa) % m2pi < (aa - ab) % m2pi:
            interval = SpatialTree.AngleInterval(aa, ab)
        else:
            interval = SpatialTree.AngleInterval(ab, aa)
    for child in spatial_tree_node.gen_children():
        interval.update(get_subtree_angle_interval(child, center))
    return interval

def rotate_subtree(spatial_tree_node, center, theta):
    """
    Modify the location values of each node in the subtree.
    @param spatial_tree_node: a spatial tree node
    @param center: the center of rotation
    @param theta: the angle (in radians) to rotate the subtree
    """
    cx, cy = center
    st = math.sin(theta)
    ct = math.cos(theta)
    for node in spatial_tree_node.gen_subtree_preorder():
        x, y = node.location
        nx = cx + (x - cx)*ct - (y - cy)*st
        ny = cy + (x - cx)*st + (y - cy)*ct
        node.location = (nx, ny)


class EqualDaylightLayoutError(Exception):
    pass


class EqualDaylightLayout:

    def __init__(self):
        self.force_straight_branches = True
        self.iterations = 5

    def do_layout(self, tree):
        """
        @param tree: something like a SpatialTree
        """
        original_root = tree.root
        if self.force_straight_branches:
            min_adjacent_nodes = 3
        else:
            min_adjacent_nodes = 2
        EqualArcLayout.do_layout(tree)
        for i in range(self.iterations):
            internal_nodes = list(tree.breadth_first())
            #random.shuffle(internal_nodes)
            for node in internal_nodes:
                adjacent_node_count = 0
                if node.parent:
                    adjacent_node_count += 1
                adjacent_node_count += len(node.children)
                if adjacent_node_count >= min_adjacent_nodes:
                    tree.reroot(node)
                    _equal_daylight(tree)
            tree.reroot(original_root)


def _equal_daylight(tree):
    """
    Rotate the child subtrees so they have equal daylight.
    """
    center = tree.root.location
    # get the intervals
    intervals = [get_subtree_angle_interval(child, center) for child in tree.root.children]
    # determine the degree to which the root node is shaded
    total_occlusion = sum(interval.get_magnitude() for interval in intervals)
    total_daylight = m2pi - total_occlusion
    if total_occlusion >= m2pi:
        raise EqualDaylightLayoutError('child subtrees occlude %f radians' % total_occlusion)
    interval_pairs = zip(intervals, intervals[1:] + [intervals[0]])
    gap_widths = [(ib.low - ia.high) % m2pi for ia, ib in interval_pairs]
    daylight_per_gap = total_daylight / float(len(tree.root.children))
    observed_cumulative_angle = 0
    expected_cumulative_angle = 0
    for child, interval, gap_width in zip(tree.root.children, intervals, gap_widths):
        angle_delta = expected_cumulative_angle - observed_cumulative_angle
        if angle_delta:
            rotate_subtree(child, center, angle_delta)
        expected_cumulative_angle += interval.get_magnitude()
        observed_cumulative_angle += interval.get_magnitude()
        expected_cumulative_angle += daylight_per_gap
        observed_cumulative_angle += gap_width


def get_test_image_format_and_string():
    import DrawTreeImage
    import Newick
    import SpatialTree
    # make a spatial tree
    tree_string = """
    ((((((((((((A:2)A:2, (B:2)B:2):3):3, (C:2.5)C:2.5):4):4, (D:3)D:3):1.5):1.5, (E:10.5)E:10.5):5):5, (((((F:2)F:2, (G:6)G:6):7):7, (H:4)H:4):6.5):6.5):6.5):6.5, (((((I:2.5)I:2.5, (J:1)J:1):15):15, (((K:5.5)K:5.5, (L:5.5)L:5.5):1):1):8.5):8.5, (M:28)M:28);"""
    tree = Newick.parse(tree_string, SpatialTree.SpatialTree)
    try:
        layout = EqualDaylightLayout()
        layout.force_straight_branches = False
        layout.do_layout(tree)
    except EqualDaylightLayoutError, e:
        print e
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
    return (image_format, image_string)


class TestEqualDaylightLayout(unittest.TestCase):

    def test(self):
        """
        Make an image in memory but do not write it to a file.
        """
        test_format, test_string = get_test_image_format_and_string()


def main():
    """
    Make an image in memory and write it to a file.
    """
    test_format, test_string = get_test_image_format_and_string()
    # write the image file
    print 'writing', filename
    filename = 'test.%s' % test_format
    fout = open(filename, 'wb')
    fout.write(test_string)
    fout.close()

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('--test', action='store_true', dest='test', default=False)
    parser.add_option('--profile', action='store_true', dest='profile', default=False, help='run some profiled unit tests')
    options, args = parser.parse_args()
    if options.test or options.profile:
        suite = unittest.TestLoader().loadTestsFromTestCase(TestEqualDaylightLayout)
        if options.profile:
            profile.run('unittest.TextTestRunner(verbosity=2).run(suite)')
        else:
            unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        main()

