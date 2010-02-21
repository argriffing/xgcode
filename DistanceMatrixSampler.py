
from StringIO import StringIO
import unittest
import math
import random

import Newick
import Clustering
import NeighborJoining
import MatrixUtil
import Util
import JC69
import NewickIO
import FelTree

def sample_infinite_sites_distances(tree, order, sequence_length):
    """
    This sampler is supposed to be somewhat fast.
    @param tree: a newick tree
    @param order: the requested order of the taxa
    @param sequence_length: the requested length of the sequences
    @return: a list of distances in the requested order
    """


def sample_infinite_alleles_sequences(tree, order, sequence_length):
    """
    This sampler is supposed to be somewhat fast.
    @param tree: a newick tree
    @param order: the requested order of the taxa
    @param sequence_length: the requested length of the sequences
    @return: a list of state sequences in the requested order
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
    # for each branch calculate the probability that no change is observed
    index_to_pequal = [-1]*n
    for i in range(1, n):
        distance = index_to_branch_length[i]
        index_to_pequal[i] = math.exp(-distance)
    # simulate nucleotide columns
    columns = []
    for iteration in range(sequence_length):
        # choose the state at the root according to the JC69 stationary distribution
        column = [0]
        # each change on a branch yields a new state without repetition
        next_state = 1
        # choose the rest of the states depending on the state of the parent
        for i in range(1, n):
            parent_state = column[index_to_parent_index[i]]
            if random.random() < index_to_pequal[i]:
                column.append(parent_state)
            else:
                column.append(next_state)
                next_state += 1
        columns.append(column)
    # convert the columns to sequences
    sequences = zip(*columns)
    # return the list of sequences in the correct order
    ordered_sequences = [sequences[name_to_index[name]] for name in order]
    return ordered_sequences

def get_infinite_alleles_ML_distance(sa, sb):
    """
    Use a closed form maximum likelihood estimator.
    @param sa: a state sequence
    @param sb: a state sequence
    @return: the estimated expected number of changes per position
    """
    assert len(sa) == len(sb)
    n = len(sa)
    differences = Util.hamming_distance(sa, sb)
    if differences == n:
        mle =  float('inf')
    else:
        mle = -math.log(float(n - differences) / float(n))
    return mle

def get_infinite_alleles_ML_distance_matrix(sequence_list):
    """
    Use a closed form maximum likelihood estimator.
    @param sequence_list: an ordered list of sequences
    @return: a row major distance matrix
    """
    return MatrixUtil.list_to_matrix(sequence_list, get_infinite_alleles_ML_distance)


class DistanceMatrixSampler:
    """
    Sample estimated distance matrices using rejection sampling.
    """

    def __init__(self, tree, ordered_names, sequence_length):
        """
        @param tree: a phylogenetic tree object
        @param ordered_names: an ordered list of names of terminal taxa
        @param sequence_length: the length of the sequences to generate (or float('inf'))
        """
        # the number of taxa should be consistent
        assert set(node.name for node in tree.gen_tips()) == set(ordered_names)
        # initialize member variables passed from the constructor
        self.tree = tree
        self.ordered_names = ordered_names
        self.sequence_length = sequence_length
        # reject branches of zero and infinite length
        self.zero_replacement = None
        self.inf_replacement = None
        # count the number of proposed and of accepted samples
        self.proposed = 0
        self.accepted = 0
        # count the number of proposed distance matrices with branch lengths of infinity and of zero
        self.proposals_with_inf = 0
        self.proposals_with_zero = 0

    def _get_computrons_per_proposal(self):
        """
        An N^k algorithm requires N^k computrons.
        @return: the number of computrons per proposal
        """
        ntaxa = len(self.ordered_names)
        if self.sequence_length == float('inf'):
            return (ntaxa ** 2)
        else:
            return (ntaxa ** 2) * self.sequence_length

    def get_completed_computrons(self):
        """
        An N^k algorithm requires N^k computrons.
        @return: the number of computrons used by the sampler so far
        """
        return self.get_completed_proposals() * self._get_computrons_per_proposal()

    def get_remaining_computrons(self, remaining_acceptances):
        """
        An N^k algorithm requires N^k computrons.
        @param remaining_acceptances: the requested number of additional samples that are accepted
        @return: an estimate of the number of remaining computrons
        """
        return self.get_remaining_proposals(remaining_acceptances) * self._get_computrons_per_proposal()

    def get_completed_proposals(self):
        """
        @return: the number of proposals that have been completed regardless of acceptance or rejection
        """
        return self.proposed

    def get_remaining_proposals(self, remaining_acceptances):
        """
        Estimate the remaining number of proposals.
        Use a prior distribution that has expectation .5 for the acceptance probability.
        @return: the estimated number of remaining iterations to accept an additional number of samples
        """
        # estimate the probability that a proposal will be accepted
        p_hat = self._get_acceptance_probability()
        # estimate the number of remaining iterations
        remaining_proposals = remaining_acceptances / p_hat
        # round the result up to an integer
        return int(remaining_proposals + 0.5)

    def _get_acceptance_probability(self):
        """
        @return: the estimated probability that a proposed sample will be accepted
        """
        # if the sequence length is infinite then we can always accept
        if self.sequence_length == float('inf'):
            return 1.0
        # use a prior acceptance probability distribution whose expected value is 0.5
        accepted_pseudocounts = 5
        proposed_pseudocounts = 10
        # calculate the posterior acceptance probability
        numerator = accepted_pseudocounts + self.accepted
        denominator = proposed_pseudocounts + self.proposed
        p_hat = float(numerator) / float(denominator)
        return p_hat

    def set_zero_replacement(self, zero_replacement):
        """
        @param zero_replacement: the value replacing a zero in the distance matrix or None to reject
        """
        self.zero_replacement = zero_replacement

    def set_inf_replacement(self, inf_replacement):
        """
        @param inf_replacement: the value replacing infinity in the distance matrix or None to reject
        """
        self.inf_replacement = inf_replacement

    def _get_proposed_sample(self):
        """
        @return: a proposed (sequence_list, distance_matrix) pair that may need to be rejected
        """
        sequence_list = JC69.sample_sequences(self.tree, self.ordered_names, self.sequence_length)
        D = JC69.get_ML_distance_matrix(sequence_list)
        return sequence_list, D

    def _process_proposed_matrix(self, D):
        """
        Modify the proposed distance matrix so that it has acceptable distances.
        This function also increments the proposals_with_inf and proposals_with_zero member variables
        @return: the modified matrix, or None if the matrix was rejected
        """
        # modify the matrix and note some of its properties
        rejected = False
        has_zero = False
        has_inf = False
        for i, row in enumerate(D):
            for j, element in enumerate(row):
                if i != j:
                    if element == 0:
                        has_zero = True
                        if self.zero_replacement is None:
                            rejected = True
                        else:
                            D[i][j] = self.zero_replacement
                    elif element == float('inf'):
                        has_inf = True
                        if self.inf_replacement is None:
                            rejected = True
                        else:
                            D[i][j] = self.inf_replacement
        if has_zero:
            self.proposals_with_zero += 1
        if has_inf:
            self.proposals_with_inf += 1
        if rejected:
            return None
        else:
            return D

    def gen_samples_or_none(self, count=None):
        """
        Yield (ordered sequence list, distance matrix) pairs or None.
        The value None is yielded if the proposal was rejected.
        @param count: the requested number of proposals or None for no bound
        """
        while True:
            # get the proposal in a form that may be yielded if it is accepted
            if self.sequence_length == float('inf'):
                # if the sequence length is infinite then no sequences are sampled
                sequence_list = None
                result = (None, self.tree.get_distance_matrix())
            else:
                result = self._get_proposed_sample()
            self.proposed += 1
            # unpack the proposal
            sequence_list, D = result
            # process the proposal, setting D to None if the proposal was rejected
            D = self._process_proposed_matrix(D)
            if D:
                self.accepted += 1
                yield result
            else:
                yield None
            # define the termination condition
            if count is not None and self.proposed >= count:
                return


class InfiniteAllelesSampler(DistanceMatrixSampler):
    """
    This class is for sampling distance matrices using the infinite alleles model.
    """

    def _get_proposed_sample(self):
        """
        @return: a proposed (None, distance_matrix) pair that may need to be rejected
        """
        sequence_list = sample_infinite_alleles_sequences(self.tree, self.ordered_names, self.sequence_length)
        D = get_infinite_alleles_ML_distance_matrix(sequence_list)
        return sequence_list, D


class TestDistanceMatrixSampler(unittest.TestCase):

    def test_finite_sequence_length(self):
        """
        Run the distance matrix sampler using a finite sequence length.
        """
        # initialize the input to the sampler
        tree_string = '(((a:0.05, b:0.05):0.15, c:0.2):0.8, x:1.0, (((m:0.05, n:0.05):0.15, p:0.2):0.8, y:1.0):1.0);'
        tree = NewickIO.parse(tree_string, FelTree.NewickTree)
        ordered_names = list(node.name for node in tree.gen_tips())
        sequence_length = 10
        # initialize the sampler, putting no limit on the number of accepted matrices or on the number of steps
        sampler = DistanceMatrixSampler(tree, ordered_names, sequence_length)
        results = list(sampler.gen_samples_or_none(100))
        self.assertEqual(len(results), 100)

    def test_modify(self):
        """
        Run the distance matrix sampler, modifying extreme distance estimates.
        """
        # initialize the input to the sampler
        tree_string = '(((a:0.05, b:0.05):0.15, c:0.2):0.8, x:1.0, (((m:0.05, n:0.05):0.15, p:0.2):0.8, y:1.0):1.0);'
        tree = NewickIO.parse(tree_string, FelTree.NewickTree)
        ordered_names = list(node.name for node in tree.gen_tips())
        sequence_length = 10
        # initialize the sampler, putting no limit on the number of accepted matrices or on the number of steps
        sampler = DistanceMatrixSampler(tree, ordered_names, sequence_length)
        # tell the sampler to modify extreme values
        sampler.set_zero_replacement('0.001')
        sampler.set_inf_replacement('20')
        # run the sampler
        results = list(sampler.gen_samples_or_none(100))

    def test_infinite_sequence_length(self):
        """
        Run the distance matrix sampler using infinite sequence length.
        """
        # initialize the input to the sampler
        tree_string = '(((a:0.05, b:0.05):0.15, c:0.2):0.8, x:1.0, (((m:0.05, n:0.05):0.15, p:0.2):0.8, y:1.0):1.0);'
        tree = NewickIO.parse(tree_string, FelTree.NewickTree)
        ordered_names = list(node.name for node in tree.gen_tips())
        sequence_length = float('inf')
        # initialize the sampler, putting no limit on the number of accepted matrices or on the number of steps
        sampler = DistanceMatrixSampler(tree, ordered_names, sequence_length)
        results = list(sampler.gen_samples_or_none(100))

    def test_runtime_estimation(self):
        """
        Run the distance matrix sampler using infinite sequence length.
        """
        # initialize the input to the sampler
        tree_string = '(((a:0.05, b:0.05):0.15, c:0.2):0.8, x:1.0, (((m:0.05, n:0.05):0.15, p:0.2):0.8, y:1.0):1.0);'
        tree = NewickIO.parse(tree_string, FelTree.NewickTree)
        ordered_names = list(node.name for node in tree.gen_tips())
        sequence_length = 10
        # initialize the sampler, putting no limit on the number of accepted matrices or on the number of steps
        sampler = DistanceMatrixSampler(tree, ordered_names, sequence_length)
        # run the sampler for a bit
        results = list(sampler.gen_samples_or_none(50))
        # get some properties of the run so far
        sampler.get_completed_computrons()
        sampler.get_remaining_computrons(100)
        sampler.get_completed_proposals()
        sampler.get_remaining_proposals(100)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestDistanceMatrixSampler)
    unittest.TextTestRunner(verbosity=2).run(suite)

