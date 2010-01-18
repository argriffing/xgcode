"""
Sample points on interlocking spirals.
This is useful for demonstrating clustering methods.
"""

import unittest
import math
import random

import numpy as np


def gen_points(npoints_per_group, noise_stddev):
    """
    Generate (x, y) pairs of points.
    @param npoints_per_group: the number of points to be sampled for each of the two labeled groups
    @param noise_stddev: the standard deviation of the noise added to the points
    """
    for x, y, label in gen_labeled_points(npoints_per_group, noise_stddev):
        yield (x, y)

def gen_labeled_points(npoints_per_group, noise_stddev):
    """
    Generate (x, y, label) triples of points, where the label is either 1 or -1
    @param npoints_per_group: the number of points to be sampled for each of the two labeled groups
    @param noise_stddev: the standard deviation of the noise added to the points
    """
    # define the number of radians in the arc that defines each sampled group
    arc_radians = .75 * 2 * math.pi
    # define the group labels
    labels = [-1, 1]
    # define the radius of the central arc that defines each group of points
    radius = 1.0
    for label in labels:
        x_center_offset = label * radius / 2
        radian_offset = (math.pi / 2) * (1 + label)
        for i in range(npoints_per_group):
            radian_sample = radian_offset + random.random() * arc_radians
            # intialize the point somewhere on a circle centered on the origin
            x = radius * math.cos(radian_sample)
            y = radius * math.sin(radian_sample)
            # move the point according to gaussian noise
            x += np.random.normal(0, noise_stddev)
            y += np.random.normal(0, noise_stddev)
            # center the point according to its group
            x += x_center_offset
            # yield the sampled point
            yield (x, y, label)


class TestSpiralSampler(unittest.TestCase):

    def test(self):
        npoints_per_group = 100
        # sample some unlabeled points
        points = list(gen_points(npoints_per_group, 1.0))
        assert len(points) == 2*npoints_per_group
        assert len(points[0]) == 2
        # sample some labeled points
        points = list(gen_labeled_points(npoints_per_group, 1.0))
        assert len(points) == 2*npoints_per_group
        assert len(points[0]) == 3


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestSpiralSampler)
    unittest.TextTestRunner(verbosity=2).run(suite)

