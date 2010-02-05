"""Explore minimum evolution Steiner trees in the sense of Edwards.

Use a numerical approximation with convergence diagnostics.
"""

# Here is what I was originally going to do:
#
#Look at trees with 4 taxa because they are the smallest trees with interesting topologies.
#First compute the embedding of the tree in 3 dimensional Euclidean space.
#Next transform the points into the (phi, r1, r2, m1, m2, h) parameterization of
#"Minimum Networks for Four Points in Space" by Rubinstein et al. 2002,
#given an {{a,b},{c,d}} split.
#Then use an iterative approach to find the two Steiner points.


from StringIO import StringIO
import math
import random

import numpy as np
import scipy

from SnippetUtil import HandlingError
import SnippetUtil
import Form
import NewickIO
import FelTree
import Euclid


g_tree_string = '((a:10, b:1):1, c:1, d:10);'



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

def get_form():
    """
    @return: a list of form objects
    """
    # define the list of form objects
    form_objects = [
            Form.MultiLine('tree_string', 'newick tree', g_tree_string)]
    return form_objects

def get_angle(x, y):
    cos_theta = np.dot(x, y) / (np.linalg.norm(x)*np.linalg.norm(y))
    return math.acos(cos_theta)

def get_random_point():
    return np.array([random.random() for i in range(3)])

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # build the newick tree from the string
    tree = NewickIO.parse(fs.tree_string, FelTree.NewickTree)
    leaf_names = list(sorted(leaf.get_name() for leaf in tree.gen_tips()))
    # assert that the newick tree has the correct set of leaves
    if leaf_names != ['a', 'b', 'c', 'd']:
        raise HandlingError('expected the tree to have leaves named {a, b, c, d}')
    # start writing the response
    out = StringIO()
    # get the distance matrix with ordered indices including all nodes in the tree
    D = np.array(tree.get_distance_matrix(leaf_names))
    # get the embedded points
    X = Euclid.edm_to_points(D)
    print >> out, 'distance matrix:'
    print >> out, D
    print >> out, 'embedded points:'
    print >> out, X
    # set up the optimization
    a, b, c, d = X.tolist()
    objective = Objective(a, b, c, d)
    s1_initial = (np.mean(X, 0) + X[0] + X[1])/3 + get_random_point()
    s2_initial = (np.mean(X, 0) + X[2] + X[3])/3 + get_random_point()
    data_initial = np.hstack([s1_initial, s2_initial])
    data_final = scipy.optimize.fmin_bfgs(objective.get_value, data_initial, fprime=objective.get_gradient, gtol=1e-10)
    s1 = data_final[:3]
    s2 = data_final[3:]
    gradient_final = objective.get_gradient(data_final)
    s1_gradient = gradient_final[:3]
    s2_gradient = gradient_final[3:]
    print >> out, 'initial random steiner point guesses:'
    print >> out, s1_initial
    print >> out, s2_initial
    print >> out, 'final steiner point estimates:'
    print >> out, s1
    print >> out, s2
    print >> out, 'each of these angles should be %f radians:' % ((2*math.pi)/3)
    print >> out, get_angle(a-s1, b-s1)
    print >> out, get_angle(b-s1, s2-s1)
    print >> out, get_angle(s2-s1, a-s1)
    print >> out, get_angle(c-s2, d-s2)
    print >> out, get_angle(d-s2, s1-s2)
    print >> out, get_angle(s1-s2, c-s2)
    print >> out, 'value of the objective function at the estimated solution:'
    print >> out, objective.get_value(data_final)
    print >> out, 'gradient of the objective function at each estimated steiner point:'
    print >> out, s1_gradient
    print >> out, s2_gradient
    # write the response
    return [('Content-Type', 'text/plain')], out.getvalue().strip()
