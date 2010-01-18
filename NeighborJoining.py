"""
Neighbor joining related algorithms.
"""

import StringIO
import unittest

import Newick
import MatrixUtil
import TreeComparison
import FelTree
import NewickIO

# This distance matrix is an example used on the Wikipedia neighbor-joining page.
g_sample_distance_matrix = [
        [ 0, 7, 11, 14],
        [ 7, 0,  6,  9],
        [11, 6,  0,  7],
        [14, 9,  7,  0]]

# This is another distance matrix used as an example.
# Primate mitochondrial DNA sequences, HindIII
# Hayasaki, Gojobori, Horai MBE 1988
# ftp.gac.edu/~mmcdermo/mcs255/j06/lectureNJ.pdf
g_mito_matrix = [
        [0, .189, .110, .113, .215],
        [.189, 0, .179, .192, .211],
        [.110, .179, 0, .094, .205],
        [.113, .192, .094, 0, .214],
        [.215, .211, .205, .214, 0]]

# These are the ordered leaf names corresponding to the mito distance matrix
g_mito_states = ['gorilla', 'orangutan', 'human', 'chimpanzee', 'gibbon']

# This is an approximation of the Q matrix for the first iteration of the mito matrix.
g_mito_matrix_q = [
        [0, -.831, -.885, -.901, -.827],
        [-.831, 0, -.822, -.808, -.983],
        [-.885, -.822, 0, -.919, -.818],
        [-.901, -.808, -.919, 0, -.816],
        [-.827, -.983, -.818, -.816, 0]]

# This is the (unrooted) tree that is supposed to be reconstructed using neighbor joining.
g_mito_tree_string = '(orangutan:0.093167, gibbon:0.11783, ((human:0.0435, chimpanzee:0.0505):0.0065, gorilla:0.058):0.0385);'

def get_Q_matrix(D):
    """
    @param D: a row major distance matrix
    @return: a row major matrix whose minimum off-diagonal defines the neighbor indices to be joined
    """
    n = len(D)
    D_star = [sum(D[i]) for i in range(n)]
    Q = []
    for i in range(n):
        row = []
        for j in range(n):
            if i == j:
                element = 0
            else:
                element = (n - 2) * D[i][j] - D_star[i] - D_star[j]
            row.append(element)
        Q.append(row)
    return Q

def get_neighbors(D):
    """
    @param D: a row major distance matrix
    @return: the pair of neighbor indices selected by the Q criterion
    """
    n = len(D)
    D_star = [sum(D[i]) for i in range(n)]
    # Use the Q matrix as a criterion to select the neighbors to join.
    Q = get_Q_matrix(D)
    q_neighbors_pairs = [(Q[i][j], (i, j)) for i in range(n) for j in range(n) if i<j]
    q, neighbors = min(q_neighbors_pairs)
    return neighbors

def do_iteration(D):
    """
    @param D: a row major distance matrix
    @return: a vector of distances to a new vertex and the pair of neighbor indices
    """
    n = len(D)
    D_star = [sum(D[i]) for i in range(n)]
    # Use the Q matrix as a criterion to select the neighbors to join.
    # The neighbor indices will be f and g.
    f, g = get_neighbors(D)
    # find the vector of distances to the new vertex that connects the neighbors
    v = [0]*n
    v[f] = D[f][g] / 2.0 + (1.0 / (2*(n-2))) * (D_star[f] - D_star[g])
    v[g] = D[g][f] / 2.0 + (1.0 / (2*(n-2))) * (D_star[g] - D_star[f])
    for i in range(n):
        if i not in (f, g):
            v[i] = (D[f][i] - v[f])/2.0 + (D[g][i] - v[g])/2.0
    return v, (f, g)

def make_tree(D, ordered_states, iteration_callback=None):
    """
    Create a newick tree from a distance matrix using neighbor joining.
    @param D: a row major distance matrix
    @param ordered_states: state names ordered according to the distance matrix
    @param iteration_callback: called with the output of each iteration call
    @return: a newick tree
    """
    # make sure that there are enough states
    if len(ordered_states) < 3:
        raise ValueError('the neighbor joining algorithm needs at least three nodes')
    # create a dictionary mapping the subtree root node serial number to a subtree
    forest = {}
    # set the current state
    index_to_serial = range(len(ordered_states))
    next_serial = len(ordered_states)
    # repeatedly pair off neighbors
    while True:
        # get the new vector of distances and the neighbor index pair
        result = do_iteration(D)
        if iteration_callback:
            # get the Q matrix for show
            Q = get_Q_matrix(D)
            # report the Q matrix and the result of the iteration
            iteration_callback(Q, result)
        v, (f, g) = result
        # create the subtree from the index pair
        root = Newick.NewickNode()
        root.serial_number = next_serial
        # determine the indices to use as branches
        if len(index_to_serial) == 3:
            branch_indices = range(3)
        else:
            branch_indices = (f, g)
        # add branches to the tree
        for index in branch_indices:
            neo = forest.pop(index_to_serial[index], None)
            if not neo:
                neo = Newick.NewickNode()
                neo.serial_number = index_to_serial[index]
            root.add_child(neo)
            neo.set_parent(root)
            neo.blen = v[index]
        # handle the terminal case
        if len(index_to_serial) == 3:
            # create the newick tree from the root node
            tree = Newick.NewickTree(root)
            # add names to the tips of the tree
            for node in tree.gen_tips():
                node.name = ordered_states[node.serial_number]
            # convert the tree to a FelTree and return it
            return NewickIO.parse(tree.get_newick_string(), FelTree.NewickTree)
        else:
            # add the subtree to the forest
            forest[next_serial] = root
            # make the next distance matrix
            next_D = []
            for i, row in enumerate(D):
                if i not in (f, g):
                    next_row = [value for j, value in enumerate(row) if j not in (f, g)]
                    next_row.append(v[i])
                    next_D.append(next_row)
            next_row = [value for j, value in enumerate(v) if j not in (f, g)]
            next_row.append(0)
            next_D.append(next_row)
            D = next_D
            # make the next serial number map
            next_index_to_serial = [value for j, value in enumerate(index_to_serial) if j not in (f, g)]
            next_index_to_serial.append(next_serial)
            index_to_serial = next_index_to_serial
            # increment the serial number
            next_serial += 1


class IterationFunctor:
    """
    This object is called for each iteration of the neighbor joining algorithm.
    """

    def __init__(self):
        self.iteration_results = []

    def __call__(self, Q, result):
        """
        @param Q: the Q matrix
        @param result: a (distance vector, index pair) pair
        """
        self.iteration_results.append((Q, result))


class TestNeighborJoining(unittest.TestCase):

    def test_get_Q_matrix(self):
        # TODO use a distance matrix with a unique best pair to simplify testing
        D = g_sample_distance_matrix
        observed_Q = get_Q_matrix(D)
        expected_Q = [
                [  0, -40, -34, -34],
                [-40,   0, -34, -34],
                [-34, -34,   0, -40],
                [-34, -34, -40,   0]]
        observed_Q_tuple = tuple(tuple(row) for row in observed_Q)
        expected_Q_tuple = tuple(tuple(row) for row in expected_Q)
        self.assertEqual(observed_Q_tuple, expected_Q_tuple)

    def test_do_iteration(self):
        D = g_mito_matrix
        f = IterationFunctor()
        make_tree(D, g_mito_states, f)
        # for the first iterations show the branch lengths
        for Q, result in f.iteration_results[:-1]:
            distance_vector, (index_a, index_b) = result
            branch_a = distance_vector[index_a]
            branch_b = distance_vector[index_b]
        # for the last iteration show the distance matrix
        Q, (distance_vector, (index_a, index_b)) = f.iteration_results[-1]

    def test_mito_matrix(self):
        D = g_mito_matrix
        n = len(D)
        observed_Q = get_Q_matrix(D)
        expected_Q = g_mito_matrix_q
        # assert that the diagonal elements of the observed Q matrix are exactly zero
        for i, row in enumerate(observed_Q):
            self.assertEqual(row[i], 0)
        # assert that the observed Q matrix is approximately equal to the expected Q matrix
        abs_tol = .001
        for i in range(n):
            for j in range(n):
                abs_delta = abs(observed_Q[i][j] - expected_Q[i][j])
                self.failUnless(abs_delta < abs_tol)
        # use neighbor joining to reconstruct the tree
        observed_tree = make_tree(D, g_mito_states)
        # load the expected tree
        expected_tree = NewickIO.parse(g_mito_tree_string, FelTree.NewickTree)
        # for the observed and expected trees calculate the induced partitions and corresponding branch lengths
        observed_partitions_and_lengths = TreeComparison.get_partitions_and_branch_lengths(observed_tree)
        expected_partitions_and_lengths = TreeComparison.get_partitions_and_branch_lengths(expected_tree)
        # the number of partitions should be the same
        self.assertEqual(len(observed_partitions_and_lengths), len(expected_partitions_and_lengths))
        # the partitions should be the same
        observed_partitions = set([part for part, length in observed_partitions_and_lengths])
        expected_partitions = set([part for part, length in expected_partitions_and_lengths])
        observed_only = observed_partitions - expected_partitions
        expected_only = expected_partitions - observed_partitions
        lines = [
                'observed partitions include: ' + str(observed_only),
                'expected partitions include: ' + str(expected_only)
                ]
        self.assertEqual(observed_partitions, expected_partitions, '\n'.join(lines))
        # corresponding partitions should have the same lengths
        observed_part_to_length = dict(observed_partitions_and_lengths)
        expected_part_to_length = dict(expected_partitions_and_lengths)
        lines = []
        for part in observed_partitions:
            observed_length = observed_part_to_length[part]
            expected_length = expected_part_to_length[part]
            abs_tol = .00001
            abs_delta = abs(observed_length - expected_length)
            if abs_delta > abs_tol:
                lines.append('partition:' + str(part))
                lines.append('observed branch length:' + str(observed_length))
                lines.append('expected branch length:' + str(expected_length))
        error_message = '\n'.join(lines)
        self.failIf(error_message, error_message)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestNeighborJoining)
    unittest.TextTestRunner(verbosity=2).run(suite)
