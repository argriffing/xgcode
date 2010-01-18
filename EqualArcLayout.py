#!/usr/bin/env python

"""
This module should only be accessed through its do_layout() function.
"""

import math

from optparse import OptionParser


def do_layout(tree):
    """
    @param tree: something like a SpatialTree
    """
    # decorate the tree with subtree_tip_count members
    for node in tree.postorder():
        if node.children:
            node.subtree_tip_count = sum(child.subtree_tip_count for child in node.children)
        else:
            node.subtree_tip_count = 1
    # create the equal arc angles
    _force_equal_arcs(tree.root, -math.pi, math.pi)
    # take off the subtree_tip_count members
    for node in tree.preorder():
        del node.subtree_tip_count
    # convert angles to coordinates
    update_locations(tree.root, (0, 0), 0)


def _force_equal_arcs(current_node, min_theta, max_theta):
    """
    Use the equal angle method to lay out the tree with non-intersecting branches.
    Define the angles in this subtree to be within the specified range.
    """
    current_node.theta = (min_theta + max_theta) / 2
    if current_node.children:
        subtree_tip_count = current_node.subtree_tip_count
        cumulative_theta = 0
        for i, child in enumerate(current_node.children):
            sub_subtree_tip_count = child.subtree_tip_count
            aliquot = (max_theta - min_theta) * sub_subtree_tip_count / float(subtree_tip_count)
            low = min_theta + cumulative_theta - current_node.theta
            cumulative_theta += aliquot
            high = min_theta + cumulative_theta - current_node.theta
            _force_equal_arcs(child, low, high)

def update_locations(current_node, last_location, last_theta):
    """
    @param last_location: an (x, y) pair or None if root
    @param last_theta: the direction we were last facing
    """
    theta = last_theta + current_node.theta
    if current_node.blen is None:
        current_node.location = last_location
    else:
        x, y = last_location
        dx = current_node.blen*math.cos(theta)
        dy = current_node.blen*math.sin(theta)
        current_node.location = (x + dx, y + dy)
    for child in current_node.children:
        update_locations(child, current_node.location, theta)


def main():
    print 'not yet...'

if __name__ == '__main__':
    parser = OptionParser()
    #parser.add_option('-v', '--verbose', action='store_true', dest='verbose', default=False)
    #parser.add_option('-o', '--output', dest='output_filename', metavar='FILE', help='output file')
    options, args = parser.parse_args()
    main()
