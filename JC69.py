"""Functions dealing with the JC69 nucleotide rate matrix.

This is the continuous time Markov model defined by Jukes and Cantor in 1969.
"""

import StringIO
import unittest
import math
import random

import MatrixUtil
import Util

def sample_xtree_sequences(root, sequence_length):
    """
    @param root: the root of a phylogenetic xtree with branch lengths
    @param sequence_length: the number of nucleotides per sequence
    @return: a list of nucleotide sequences ordered by the xtree leaf labels
    """
    ntips = len(root.get_labeled_vertices())
    sequences = [None]*ntips
    id_to_sequence = {}
    for vertex in root.get_preorder_vertices():
        if vertex.branch:
            letters = []
            p_randomize = 1 - math.exp(-(4.0/3.0)*vertex.branch.length)
            for parent_nt in id_to_sequence[id(vertex.branch.parent)]:
                if random.random() < p_randomize:
                    letters.append(random.choice('ACGT'))
                else:
                    letters.append(parent_nt)
        else:
            letters = [random.choice('ACGT') for i in range(sequence_length)]
        sequence = ''.join(letters)
        if vertex.has_label():
            sequences[vertex.label] = sequence
        else:
            id_to_sequence[id(vertex)] = sequence
    return sequences

def sample_sequences(tree, order, sequence_length):
    """
    This sampler is supposed to be somewhat fast.
    @param tree: a newick tree
    @param order: the requested order of the taxa
    @param sequence_length: the requested length of the sequences
    @return: a list of nucleotide sequences in the requested order
    """
    # assert that the taxa named in the order are the same as the names of the tips
    sorted_tree_tip_names = tuple(sorted(node.get_name() for node in tree.gen_tips()))
    sorted_order_names = tuple(sorted(order))
    assert sorted_tree_tip_names == sorted_order_names
    # get the number of nodes in the tree, including internal nodes
    n = len(list(tree.preorder()))
    # make a preorder array representation of the tree, setting the defaults for the root
    index_to_parent_index = [-1]*n
    index_to_branch_length = [-1]*n
    # map the node id and the node name to the index
    name_to_index = {}
    id_to_index = {}
    for i, node in enumerate(tree.preorder()):
        name_to_index[node.get_name()] = i
        id_to_index[id(node)] = i
        if i:
            index_to_parent_index[i] = id_to_index[id(node.get_parent())]
            index_to_branch_length[i] = node.get_branch_length()
    # for each branch calculate the probability that four way randomization is not needed
    index_to_pequal = [-1]*n
    for i in range(1, n):
        distance = index_to_branch_length[i]
        index_to_pequal[i] = math.exp(-(4.0/3.0)*distance)
    # simulate nucleotide columns
    columns = []
    for iteration in range(sequence_length):
        # choose the state at the root according to the JC69 stationary distribution
        column = [random.choice('ACGT')]
        # choose the rest of the states depending on the state of the parent
        for i in range(1, n):
            parent_state = column[index_to_parent_index[i]]
            if random.random() < index_to_pequal[i]:
                column.append(parent_state)
            else:
                column.append(random.choice('ACGT'))
        columns.append(column)
    # convert the columns to sequences
    sequences = zip(*columns)
    # return the list of sequences in the correct order
    ordered_sequences = [''.join(sequences[name_to_index[name]]) for name in order]
    return ordered_sequences

def get_ML_distance(sa, sb):
    """
    Use a closed form maximum likelihood estimator.
    @param sa: a nucleotide sequence
    @param sb: a nucleotide sequence
    @return: the estimated expected number of changes per position
    """
    assert len(sa) == len(sb)
    assert set(sa+sb) <= set('ACGT')
    n = len(sa)
    n2 = Util.hamming_distance(sa, sb)
    n1 = n - n2
    if n2 >= 3 * n1:
        return float('inf')
    mle = -(3.0/4.0) * math.log((3.0*n1 - n2) / (3.0*n1 + 3.0*n2))
    return mle

def get_ML_distance_matrix(sequence_list):
    """
    Use a closed form maximum likelihood estimator.
    @param sequence_list: an ordered list of nucleotide sequences
    @return: a row major distance matrix
    """
    return MatrixUtil.list_to_matrix(sequence_list, get_ML_distance)


class TestJC69(unittest.TestCase):
    
    def test_placeholder(self):
        pass


if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('--test', action='store_true', dest='test', default=False)
    options, args = parser.parse_args()
    if options.test:
        suite = unittest.TestLoader().loadTestsFromTestCase(TestJC69)
        unittest.TextTestRunner(verbosity=2).run(suite)

