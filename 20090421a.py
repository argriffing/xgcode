"""Examine the relationship between an eigenvalue and the quality of its corresponding spectral split.

Phylogenetic trees are sampled by random agglomeration of some number of taxa.
Distance matrices are sampled from these trees.
Trees are estimated from these sampled distance matrices.
R-readable output suitable for making boxplots is generated,
with lines corresponding to each combination of three normalizations and two split outcomes.
The two outcomes are valid and invalid splits.
The three normalizations are: unnormalized, normalized by the sum of the first two eigenvalues,
and normalized by the sum of all eigenvalues.
"""

import StringIO
import time
import math
import random
import optparse

import numpy
from numpy import linalg

from SnippetUtil import HandlingError
import MatrixUtil
import TreeSampler
import BuildTreeTopology
import BranchLengthSampler
import JC69
import Euclid
import Form

class InfiniteDistanceError(Exception): pass
class ZeroDistanceError(Exception): pass

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.Integer('length', 'sequence length for distance matrix sampling', 1000, low=10, high=20000),
            Form.Integer('ntaxa', 'number of taxa for distance matrix sampling', 20, low=4, high=20),
            Form.RadioGroup('tree_sampling', 'branch length distribution', [
                Form.RadioItem('pachter_length', str(BranchLengthSampler.Pachter()), True),
                Form.RadioItem('exponential_length', str(BranchLengthSampler.Exponential())),
                Form.RadioItem('uniform_length_a', str(BranchLengthSampler.UniformA())),
                Form.RadioItem('uniform_length_b', str(BranchLengthSampler.UniformB()))])]
    return form_objects

def sample_distance_matrix(xtree_root, sequence_length):
    sequences = JC69.sample_xtree_sequences(xtree_root, sequence_length)
    nsequences = len(sequences)
    pairwise_mismatch_count = numpy.zeros((nsequences, nsequences))
    for i, sa in enumerate(sequences):
        for j, sb in enumerate(sequences):
            if i < j:
                nmismatches = sum(1 for a, b in zip(sa, sb) if a != b)
                if not nmismatches:
                    raise ZeroDistanceError()
                if nmismatches * 4 >= sequence_length * 3:
                    raise InfiniteDistanceError()
                pairwise_mismatch_count[i][j] = nmismatches
    D = numpy.zeros_like(pairwise_mismatch_count)
    for i in range(nsequences):
        for j in range(nsequences):
            if i < j:
                d_raw = pairwise_mismatch_count[i][j] / float(sequence_length)
                b = 0.75
                d_mle = -b*math.log(1 - d_raw/b)
                D[i][j] = d_mle
                D[j][i] = d_mle
    return D


def list_to_c(arr):
    """
    @param arr: an array of real numbers
    @return: an R string
    """
    return 'c(' + ', '.join(str(x) for x in arr) + ')'


class Builder:

    def __init__(self):
        self.true_splits = None
        self.eigenvalues = None
        self.succeeded_unnormalized = []
        self.succeeded_normalized_2 = []
        self.succeeded_normalized_n = []
        self.failed_unnormalized = []
        self.failed_normalized_2 = []
        self.failed_normalized_n = []

    def get_unnormalized(self):
        """
        @return: the unnormalized greatest eigenvalue
        """
        return max(self.eigenvalues)
    
    def get_normalized_2(self):
        """
        Normalization is by the sum of the greatest two eigenvalues.
        @return: the normalized greatest eigenvalue
        """
        top_pair = list(sorted(self.eigenvalues))[-2:]
        return max(self.eigenvalues) / sum(top_pair)

    def get_normalized_n(self):
        """
        Normalization is by the sum of all eigenvalues.
        @return: the normalized greatest eigenvalue
        """
        return max(self.eigenvalues) / sum(self.eigenvalues)

    def evaluate(self, true_splits, D_estimated):
        """
        @param true_splits: the set of all full splits implied by the true tree
        @param D_estimated: the estimated distance matrix
        """
        self.true_splits = true_splits
        BuildTreeTopology.get_splits(D_estimated, self.split_function, BuildTreeTopology.update_using_laplacian, self.on_label_split)

    def split_function(self, D):
        """
        Split the distance matrix using signs of an eigenvector of -HDH/2.
        If a degenerate split is found then a DegenerateSplitException is raised.
        @param D: the distance matrix
        @return: a set of two index sets defining a split of the indices
        """
        try:
            # get the matrix whose eigendecomposition is of interest
            HSH = Euclid.edm_to_dccov(D)
            # get the eigendecomposition
            eigenvalues, V_T = linalg.eigh(HSH)
            eigenvectors = V_T.T.tolist()
            # save the eigenvalues for reporting
            self.eigenvalues = eigenvalues
            # get the eigenvector of interest
            w, v = max(zip(eigenvalues, eigenvectors))
            # get the indices with positive eigenvector valuations
            n = len(D)
            positive = frozenset(i for i, x in enumerate(v) if x > 0)
            nonpositive = frozenset(set(range(n)) - positive)
            # check for a degenerate split
            for index_set in (positive, nonpositive):
                assert len(index_set) > 0
            for index_set in (positive, nonpositive):
                if len(index_set) == 1:
                    index, = index_set
                    raise BuildTreeTopology.DegenerateSplitException(index)
            return frozenset((positive, nonpositive))
        except BuildTreeTopology.DegenerateSplitException, e:
            self.eigenvalues = None
            return BuildTreeTopology.split_nj(D)

    def on_label_split(self, label_split):
        """
        This is a callback function that is called when a label split is calculated.
        At this point the true splits of the tree are known.
        The eigenvalues associated with the distance matrix split that generated this
        label split are also known.
        @param label_split: a reported split of the labels
        """
        # If the eigenvalues are not available then the last split
        # fell back to neighbor joining.
        if self.eigenvalues is not None:
            if label_split in self.true_splits:
                self.succeeded_unnormalized.append(self.get_unnormalized())
                self.succeeded_normalized_2.append(self.get_normalized_2())
                self.succeeded_normalized_n.append(self.get_normalized_n())
            else:
                self.failed_unnormalized.append(self.get_unnormalized())
                self.failed_normalized_2.append(self.get_normalized_2())
                self.failed_normalized_n.append(self.get_normalized_n())

    def get_R_strings(self):
        results = [
                'succeeded.unnormalized <- ' + list_to_c(self.succeeded_unnormalized),
                'succeeded.normalized_2 <- ' + list_to_c(self.succeeded_normalized_2),
                'succeeded.normalized_n <- ' + list_to_c(self.succeeded_normalized_n),
                'failed.unnormalized <- ' + list_to_c(self.failed_unnormalized),
                'failed.normalized_2 <- ' + list_to_c(self.failed_normalized_2),
                'failed.normalized_n <- ' + list_to_c(self.failed_normalized_n)]
        return results


def process(ntaxa, length, nseconds, branch_length_sampler):
    """
    @param ntaxa: the number of taxa in the sampled trees
    @param length: the length of sequences used to sample the distance matrix
    @param nseconds: allow this many seconds to run or None to run forever
    @param branch_length_sampler: a functor that returns a branch length and has a string cast
    @return: a multi-line string that summarizes the results
    """
    start_time = time.time()
    # initialize the builder object
    builder = Builder()
    # track the number of samples that failed for various reasons
    n_zero_errors = 0
    n_infinite_errors = 0
    n_failed_spectral_splits = 0
    # do a bunch of reconstructions of sampled distance matrices
    try:
        while True:
            elapsed_time = time.time() - start_time
            if nseconds and elapsed_time > nseconds:
                break
            # sample the tree topology and get its set of implied full label splits
            tree = TreeSampler.sample_agglomerated_tree(ntaxa)
            true_splits = tree.get_nontrivial_splits()
            # sample the branch lengths
            for branch in tree.get_branches():
                branch.length = branch_length_sampler()
            try:
                D = sample_distance_matrix(tree, length)
                # record information about the splits
                builder.evaluate(true_splits, D)
            except InfiniteDistanceError, e:
                n_infinite_errors += 1
            except ZeroDistanceError, e:
                n_zero_errors += 1
            except BuildTreeTopology.InvalidSpectralSplitException, e:
                n_failed_spectral_splits += 1
    except KeyboardInterrupt, e:
        pass
    # make the response
    out = StringIO.StringIO()
    print >> out, '#', elapsed_time, 'seconds of run time'
    print >> out, '#', ntaxa, 'taxa per tree'
    print >> out, '#', branch_length_sampler
    print >> out, '#', length, 'nucleotides per sequence'
    print >> out, '#', n_zero_errors, 'samples were rejected because of an estimated distance of zero'
    print >> out, '#', n_infinite_errors, 'samples were rejected because of an estimated distance of infinity'
    print >> out, '#', n_failed_spectral_splits, 'distance matrices induced an intermediate distance matrix with no spectral split'
    for s in builder.get_R_strings():
        print >> out, s
    return out.getvalue().strip()

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    nseconds = 2
    length = fs.length
    ntaxa = fs.ntaxa
    # define the branch length sampler
    if fs.pachter_length:
        branch_length_sampler = BranchLengthSampler.Pachter()
    elif fs.exponential_length:
        branch_length_sampler = BranchLengthSampler.Exponential()
    elif fs.uniform_length_a:
        branch_length_sampler = BranchLengthSampler.UniformA()
    elif fs.uniform_length_b:
        branch_length_sampler = BranchLengthSampler.UniformB()
    response_text = process(ntaxa, length, nseconds, branch_length_sampler)
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, response_text

def main(options):
    assert 0 <= options.nseconds
    assert 4 <= options.ntaxa <= 20
    assert 10 <= options.length <= 10000
    branch_length_sampler = BranchLengthSampler.Pachter()
    print process(options.ntaxa, options.length, options.nseconds, branch_length_sampler)

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('--nseconds', dest='nseconds', type='int', default=0, help='seconds to run or 0 to run until ctrl-c')
    parser.add_option('--ntaxa', dest='ntaxa', type='int', default=20, help='number of taxa in each sampled tree topology')
    parser.add_option('--length', dest='length', type='int', default=200, help='nucleotide sequence length for distance matrix sampling')
    options, args = parser.parse_args()
    main(options)
