"""Given a newick tree, find its center according to some criterion.
"""

from StringIO import StringIO
import math
import random

import numpy as np

from SnippetUtil import HandlingError
import MatrixUtil
import NewickIO
import FelTree
import MST
import Form
import FormOut

#FIXME this snippet duplicates library code

def get_form():
    """
    @return: the body of a form
    """
    # define the default tree string
    default_tree_string = '(a:1, (b:2, d:5):1, c:4);'
    # define the form objects
    return [Form.MultiLine('tree', 'newick tree', default_tree_string)]

def get_form_out():
    return FormOut.Report()

def get_laplacian_pseudo_inverse(distance_matrix):
    """
    @param distance_matrix: a row major distance matrix
    @return: the pseudo inverse laplacian matrix as a numpy array
    """
    n = len(distance_matrix)
    M = np.array(distance_matrix)
    P = np.eye(n) - np.ones((n,n))/n
    L_pinv = - 0.5 * np.dot(P, np.dot(M, P))
    return L_pinv

def get_eigendecomposition(M):
    """
    @param M: a numpy array
    @return: the eigenvalues and the eigenvectors
    """
    w, v = np.linalg.eigh(M)
    eigenvalues = w
    eigenvectors = v.T
    return eigenvalues, eigenvectors

def gen_euclidean_points_from_eigendecomposition(w, v):
    """
    Yields euclidean points.
    @param w: eigenvalues
    @param v: eigenvectors
    """
    n = len(w)
    for i in range(n):
        z = [v[j][i] * math.sqrt(w[j]) for j in range(n)]
        yield z

def gen_euclidean_points(distance_matrix):
    """
    Convert an N by N distance matrix into a list of N N-dimensional vectors.
    Yields euclidean points that should give the same distance matrix
    if the original distance were among taxa.
    The input distance matrix is a row major matrix of distances
    between taxa on a tree.
    @param distance_matrix: distances between tree taxa
    """
    # get the pseudo inverse laplacian matrix
    L_pinv = get_laplacian_pseudo_inverse(distance_matrix)
    # get the eigendecomposition of the pseudo inverse laplacian matrix
    w, v = get_eigendecomposition(L_pinv)
    """
    print >> out, 'eigenvalues of the pseudo inverse of the laplacian:'
    print >> out, eigenvalues
    print >> out, 'eigenvectors of the pseudo inverse of the laplacian:'
    print >> out, eigenvectors
    """
    for point in gen_euclidean_points_from_eigendecomposition(w, v):
        yield point

def get_euclidean_distance_matrix(points):
    """
    @param points: a sequence of euclidean points
    @return: a distance matrix
    """
    D = []
    for point_a in points:
        row = []
        for point_b in points:
            distance = math.sqrt(sum((b - a)**2
                for a, b in zip(point_a, point_b)))
            row.append(distance)
        D.append(row)
    return D

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # get the tree
    tree = NewickIO.parse(fs.tree, FelTree.NewickTree)
    states = list(sorted(node.name for node in tree.gen_tips()))
    n = len(states)
    # start to prepare the reponse
    out = StringIO()
    # get the distance matrix
    distance_matrix = tree.get_distance_matrix(states)
    # get the equivalent euclidean points
    z_points = list(gen_euclidean_points(distance_matrix))
    # get the centroid
    centroid = [sum(values)/n for values in zip(*z_points)]
    # get the resistance distances between the centroid and each point
    #volume = -sum(L[i][j] for i in range(n) for j in range(n) if i != j)
    #volume *= (4.0 / 4.3185840708)
    #volume = 1
    """
    print >> out, 'distances to the first point:'
    for z in z_points:
        print >> out, sum((a-b)**2 for a, b in zip(z, z_points[0]))
    print >> out, 'distances to the centroid:'
    for z in z_points:
        print >> out, sum((a-b)**2 for a, b in zip(z, centroid))
    """
    print >> out, 'distances to the virtual center of the tree:'
    origin = [0 for i in range(n)]
    for z in z_points:
        print >> out, sum((a-b)**2 for a, b in zip(z, origin))
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()

def hard_coded_analysis_a():
    tree_string = '(a:1, (b:2, d:5):1, c:4);'
    tree = NewickIO.parse(tree_string, FelTree.NewickTree)
    states = []
    id_list = []
    for state, id_ in sorted((node.name, id(node))
            for node in tree.gen_tips()):
        id_list.append(id_)
        states.append(state)
    for node in tree.gen_internal_nodes():
        id_list.append(id(node))
        states.append('')
    n = len(states)
    for method in ('tips', 'full'):
        # get the distance matrix from the tree
        if method == 'tips':
            print 'leaves only:'
            distance_matrix = tree.get_distance_matrix(states)
        else:
            print 'leaves and internal nodes:'
            distance_matrix = tree.get_full_distance_matrix(id_list)
        print 'distance matrix from the tree:'
        print MatrixUtil.m_to_string(distance_matrix)
        # get the equivalent euclidean points
        z_points = list(gen_euclidean_points(distance_matrix))
        for state, point in zip(states, z_points):
            print state, point
        # get the distance matrix from the transformed points
        print 'distance matrix from the transformed points:'
        distance_matrix = get_euclidean_distance_matrix(z_points)
        print MatrixUtil.m_to_string(distance_matrix)
        print


class MyObjective:

    def __init__(self, z_points):
        """
        @param z_points: N N-dimensional vectors, one for each taxon
        """
        self.z_points = z_points
        self.best = None

    def __call__(self, v):
        """
        @param v: a vector
        @return: the cost of this vector
        """
        n = len(self.z_points)
        assert len(v) == 2*n
        # first split v into its two constituent vectors
        va, vb = v[:n], v[n:]
        # get the augmented list of points
        augmented_points = self.z_points + [va, vb]
        # get the distance matrix
        distance_matrix = get_euclidean_distance_matrix(augmented_points)
        # define the vertices of the graph
        V = range(len(augmented_points))
        # define the weighted edges of the graph
        E = []
        for i in range(len(distance_matrix)):
            for j in range(len(distance_matrix)):
                if i < j:
                    distance = distance_matrix[i][j]
                    weight = distance ** 2
                    edge = (weight, i, j)
                    E.append(edge)
        # find the minimum spanning tree
        T = MST.kruskal(V, E)
        # find the total weight of the minimum spanning tree
        total_weight = sum(weight for weight, a, b in T)
        # track the best found so far
        state = (total_weight, T)
        if self.best is None:
            self.best = state
        self.best = min(self.best, state)
        # return the total weight that we want minimized
        return total_weight


def hard_coded_analysis_b():
    """
    Numerically search for the power 2 steiner points.
    """
    # make a distance matrix where the order is alphabetical with the states
    tree_string = '(a:1, (b:2, d:5):1, c:4);'
    tree = NewickIO.parse(tree_string, FelTree.NewickTree)
    states = list(sorted(node.name for node in tree.gen_tips()))
    distance_matrix = tree.get_distance_matrix(states)
    # get the pseudo inverse laplacian matrix
    L_pinv = get_laplacian_pseudo_inverse(distance_matrix)
    # get the eigendecomposition of the pseudo inverse laplacian matrix
    eigenvalues, eigenvectors = get_eigendecomposition(L_pinv)
    print 'eigenvalues of the pseudo inverse of the laplacian:'
    print eigenvalues
    # each taxon gets a transformed point
    z_points = list(gen_euclidean_points_from_eigendecomposition(
        eigenvalues, eigenvectors))
    # initialize the objective function
    objective = MyObjective(z_points)
    # initialize a couple of steiner points
    n = len(states)
    va = [random.random() for i in range(n)]
    vb = [random.random() for i in range(n)]
    # define the initial guess
    x0 = va + vb
    # do the optimization
    result = optimize.fmin(objective, x0)
    print result
    print objective.best



def main():
    hard_coded_analysis_b()

if __name__ == '__main__':
    main()
