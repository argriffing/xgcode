"""Analyze a tetrahedron in a paper by Cavalli-Sforza and A.W.F. Edwards.

Analyze the tetrahedron from the paper 
"Phylogenetic Analysis Models and Estimation Procedures."
tetrahedron.
In the paper,
Cavalli-Sforza and A.W.F Edwards give branch lengths for three trees
that have been derived from a tetrahedron whose four vertices are
{(K)orean, E(S)kimo, (B)antu, and (E)nglish}.
The topologies of the estimated trees are each ((K, S), (B, E)).
The three estimation methods are
{"additive tree", "minimum evolution", "maximum likelihood"}.
For the additive tree,
the internal vertices extend into dimensions other than the three
in which the given tetrahedron resides.
The other two trees reside entirely in the 3D space of the tetrahedron.
"""


from StringIO import StringIO
import math

import numpy as np
import scipy

from SnippetUtil import HandlingError
import SnippetUtil
import Form
import FormOut
import NewickIO
import NeighborJoining
import FelTree
import Steiner


def get_form():
    """
    @return: a list of form objects
    """
    # define the list of form objects
    form_objects = [
            Form.Matrix('tetrahedron',
                'tetrahedron with an underlying ((1,2),(3,4)) topology',
                Steiner.g_X)]
    return form_objects

def get_form_out():
    return FormOut.Report()

def do_steiner_analysis(X):
    a, b, c, d = X.tolist()
    objective = Steiner.Objective(a, b, c, d)
    s1_initial = (np.mean(X, 0) + X[0] + X[1])/3 + Steiner.get_random_point()/10
    s2_initial = (np.mean(X, 0) + X[2] + X[3])/3 + Steiner.get_random_point()/10
    data_initial = np.hstack([s1_initial, s2_initial])
    data_final = scipy.optimize.fmin_bfgs(objective.get_value, data_initial, fprime=objective.get_gradient, gtol=1e-10)
    s1 = data_final[:3]
    s2 = data_final[3:]
    gradient_final = objective.get_gradient(data_final)
    s1_gradient = gradient_final[:3]
    s2_gradient = gradient_final[3:]
    out = StringIO()
    print >> out, 'initial random steiner point guesses:'
    print >> out, s1_initial
    print >> out, s2_initial
    print >> out, 'final steiner point estimates:'
    print >> out, s1
    print >> out, s2
    print >> out, 'each of these angles should be %f radians:' % ((2*math.pi)/3)
    print >> out, Steiner.get_angle(a-s1, b-s1)
    print >> out, Steiner.get_angle(b-s1, s2-s1)
    print >> out, Steiner.get_angle(s2-s1, a-s1)
    print >> out, Steiner.get_angle(c-s2, d-s2)
    print >> out, Steiner.get_angle(d-s2, s1-s2)
    print >> out, Steiner.get_angle(s1-s2, c-s2)
    print >> out, 'value of the objective function at the estimated solution:'
    print >> out, objective.get_value(data_final)
    print >> out, 'gradient of the objective function at each estimated steiner point:'
    print >> out, s1_gradient
    print >> out, s2_gradient
    print >> out, 'branch lengths:'
    print >> out, 'branch 0:', np.linalg.norm(a-s1)
    print >> out, 'branch 1:', np.linalg.norm(b-s1)
    print >> out, 'branch 2:', np.linalg.norm(c-s2)
    print >> out, 'branch 3:', np.linalg.norm(d-s2)
    print >> out, 'central branch:', np.linalg.norm(s2-s1)
    return out.getvalue().strip()

def do_distance_analysis(X):
    # get the matrix of squared distances
    labels = list('0123')
    # reconstruct the matrix of Euclidean distances from a tree
    D_sqrt = np.array([[np.linalg.norm(y-x) for x in X] for y in X])
    sqrt_tree = NeighborJoining.make_tree(D_sqrt, labels)
    sqrt_tree_string = NewickIO.get_newick_string(sqrt_tree)
    sqrt_feltree = NewickIO.parse(sqrt_tree_string, FelTree.NewickTree)
    D_sqrt_reconstructed = np.array(sqrt_feltree.get_distance_matrix(labels))
    # reconstruct the matrix of squared Euclidean distances from a tree
    D = D_sqrt**2
    tree = NeighborJoining.make_tree(D, labels)
    tree_string = NewickIO.get_newick_string(tree)
    feltree = NewickIO.parse(tree_string, FelTree.NewickTree)
    D_reconstructed = np.array(feltree.get_distance_matrix(labels))
    # start writing
    out = StringIO()
    # matrix of Euclidean distances and its reconstruction from a tree
    print >> out, 'matrix of Euclidean distances between tetrahedron vertices:'
    print >> out, D_sqrt
    print >> out, 'neighbor joining tree constructed from D = non-squared Euclidean distances (unusual):'
    print >> out, sqrt_tree_string
    print >> out, 'distance matrix implied by this tree:'
    print >> out, D_sqrt_reconstructed
    # matrix of squared Euclidean distances and its reconstruction from a tree
    print >> out, 'matrix of squared distances between tetrahedron vertices:'
    print >> out, D
    print >> out, 'neighbor joining tree constructed from D = squared Euclidean distances (normal):'
    print >> out, tree_string
    print >> out, 'distance matrix implied by this tree:'
    print >> out, D_reconstructed
    return out.getvalue().strip()

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    X = fs.tetrahedron
    # start writing the response
    out = StringIO()
    print >> out, do_steiner_analysis(X)
    print >> out, do_distance_analysis(X)
    # write the response
    return [('Content-Type', 'text/plain')], out.getvalue().strip()
