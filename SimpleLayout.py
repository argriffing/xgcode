#!/usr/bin/env python

"""
This is an implementation of the equal angle layout algorithm.
It is described on pages 578-580 of Inferring Phylogenies.
This module should only be accessed through its do_layout() function.
"""

import math

from optparse import OptionParser


def do_layout(tree):
    """
    @param tree: something like a SpatialTree
    """
    _force_equal_angles(tree.root, -math.pi, math.pi)
    update_locations(tree.root, (0, 0), 0)

def _force_equal_angles(current_node, min_theta, max_theta):
    """
    Use the equal angle method to lay out the tree with non-intersecting branches.
    Define the angles in this subtree to be within the specified range.
    """
    current_node.theta = (min_theta + max_theta) / 2
    if current_node.children:
        aliquot = (max_theta - min_theta) / len(current_node.children)
        for i, child in enumerate(current_node.children):
            low = min_theta + i*aliquot - current_node.theta
            high = min_theta + (i+1)*aliquot - current_node.theta
            _force_equal_angles(child, low, high)

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
