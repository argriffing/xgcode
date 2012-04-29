"""
Make a minimum Steiner tree.
"""


import unittest
import math
import random

import numpy as np
import scipy
from scipy import optimize


# Table 3 in Cavalli-Sforza and Edwards (1967),
# reorded as Eskimo, Korean, Bantu, English as in table 7.
g_X = np.array([
    [0, 0, 0],
    [.0379, .0428, .4885],
    [.9136, 0, 0],
    [.4695, .3895, 0],
    ])


class Objective:

    def __init__(self, a, b, c, d):
        """
        @param a: a tip on the left side of the steiner topology
        @param b: a tip on the left side of the steiner topology
        @param c: a tip on the right side of the steiner topology
        @param d: a tip on the right side of the steiner topology
        """
        self.a = a
        self.b = b
        self.c = c
        self.d = d

    def get_value(self, data):
        """
        The idea is to minimize the value returned by this function.
        @param data: a vector representing a point in the 6d search space
        @return: the value of the objective function
        """
        s1 = data[:3]
        s2 = data[3:]
        return get_tree_length(self.a, self.b, self.c, self.d, s1, s2)

    def get_gradient(self, data):
        """
        Compute the s1 gradient and s2 gradient separately and concatenate them.
        @param data: a vector representing a point in the 6d search space
        @return: the gradient of the objective function
        """
        a = self.a
        b = self.b
        c = self.c
        d = self.d
        s1 = data[:3]
        s2 = data[3:]
        s1_gradient = normalized(s1 - a) + normalized(s1 - s2) + normalized(s1 - b)
        s2_gradient = normalized(s2 - c) + normalized(s2 - s1) + normalized(s2 - d)
        return np.hstack([s1_gradient, s2_gradient])


def normalized(v):
    return v / np.linalg.norm(v)

def get_tree_length(a, b, c, d, s1, s2):
    """
    @param a: a tip on the left side of the steiner topology
    @param b: a tip on the left side of the steiner topology
    @param c: a tip on the right side of the steiner topology
    @param d: a tip on the right side of the steiner topology
    @param s1: a steiner point on the left side of the steiner topology
    @param s2: a steiner point on the right side of the steiner topology
    """
    edges = [(a, s1), (b, s1), (s1, s2), (s2, c), (s2, d)]
    return sum(np.linalg.norm(y-x) for x, y in edges)

def get_angle(x, y):
    cos_theta = np.dot(x, y) / (np.linalg.norm(x)*np.linalg.norm(y))
    return math.acos(cos_theta)

def get_random_point():
    return np.array([random.random() for i in range(3)])


class TestSteiner(unittest.TestCase):

    def test_a(self):
        # set up the embedded points and the implicity topology of the steiner tree
        X = g_X
        a, b, c, d = X.tolist()
        objective = Objective(a, b, c, d)
        # get the estimated Steiner points and the associated gradient which should be zero
        s1_initial = (np.mean(X, 0) + X[0] + X[1])/3 + get_random_point()/10
        s2_initial = (np.mean(X, 0) + X[2] + X[3])/3 + get_random_point()/10
        data_initial = np.hstack([s1_initial, s2_initial])
        data_final = scipy.optimize.fmin_bfgs(objective.get_value, data_initial, fprime=objective.get_gradient, gtol=1e-10, disp=False)
        s1 = data_final[:3]
        s2 = data_final[3:]
        gradient_final = objective.get_gradient(data_final)
        s1_gradient = gradient_final[:3]
        s2_gradient = gradient_final[3:]
        # compare the Steiner points to expected values
        s1_expected = np.array([0.12205791, 0.07466552, 0.13937451])
        s2_expected = np.array([0.47098798, 0.26563362, 0.02757105])
        self.assertTrue(np.allclose(s1, s1_expected))
        self.assertTrue(np.allclose(s2, s2_expected))
        # compare the gradients to zero
        self.assertTrue(np.allclose(s1_gradient, np.zeros(3)))
        self.assertTrue(np.allclose(s2_gradient, np.zeros(3)))


if __name__ == '__main__':
    unittest.main()
