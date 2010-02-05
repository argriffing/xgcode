"""Compare the quality of reconstructed trees using different methods.

Phylogenetic trees are sampled by random agglomeration of some number of taxa.
Distance matrices are sampled from these trees.
Trees are estimated from these sampled distance matrices.
Each row of data output corresponds to
a sampled tree,
a nucleotide alignment sampled on the tree using the Jukes-Cantor model,
a matrix of Jukes-Cantor corrected distances,
and the trees that have been estimated from this distance matrix.
Each distance matrix is checked for the Atteson condition with respect to the true tree.
For each method,
the Robinson-Foulds distance
between the estimated tree and the true tree is calculated.
"""

from StringIO import StringIO
import time
import math
import random
import optparse

import numpy

from SnippetUtil import HandlingError
import MatrixUtil
import BranchLengthSampler
import TreeSampler
import BuildTreeTopology
import Xtree
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
                Form.RadioItem('uniform_length_b', str(BranchLengthSampler.UniformB()))]),
            Form.CheckGroup('method_options', 'compare these methods', [
                Form.CheckItem('nj', 'neighbor joining', True),
                Form.CheckItem('modified_nj', 'neighbor joining with laplacian updates', True),
                Form.CheckItem('all_spectral', 'spectral splitting with laplacian updates and neighbor joining fallback', True),
                Form.CheckItem('one_spectral', 'an initial spectral split followed by neighbor joining', True)]),
            Form.RadioGroup('delivery', 'delivery', [
                Form.RadioItem('inline', 'view as text', True),
                Form.RadioItem('attachment', 'download as an R table')])]
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


class ScatterPoint:
    """
    Each instance corresponds to a point on the scatter plot.
    """

    def __init__(self, atteson, nj_error, modified_nj_error, all_spectral_error, one_spectral_error):
        """
        Save the Robinson-Foulds errors.
        @param atteson: True iff the distance matrix is Atteson with respect to the original tree
        @param nj_error: the RF distance from the true tree to the neighbor joining tree
        @param modified_nj_error: the distance to the tree estimated using nj splitting and Laplacian updating
        @param all_spectral_error: the distance to the all-spectral-cut (with fallback) tree
        @param one_spectral_error: the distance to the single-spectral-cut-then-nj tree
        """
        self.atteson = atteson
        self.nj_error = nj_error
        self.modified_nj_error = modified_nj_error
        self.all_spectral_error = all_spectral_error
        self.one_spectral_error = one_spectral_error


class SplitFunctor:
    """
    This split functor splits large matrices differently than small matrices.
    """

    def __init__(self, large_matrix_size):
        """
        @param large_matrix_size: distance matrices at least this big will be split spectrally
        """
        self.large_matrix_size = large_matrix_size

    def __call__(self, D):
        """
        @param D: the distance matrix
        @return: a set of two index sets defining a split of the indices
        """
        if len(D) < self.large_matrix_size:
            return BuildTreeTopology.split_nj(D)
        else:
            return BuildTreeTopology.split_using_eigenvector_with_nj_fallback(D)


class UpdateFunctor:
    """
    This update functor updates large matrices differently than small matrices.
    """

    def __init__(self, large_matrix_size):
        """
        @param large_matrix_size: distance matrices at least this big will be updated using the laplacian
        """
        self.large_matrix_size = large_matrix_size

    def __call__(self, D, index_set):
        """
        @param D: the distance matrix
        @param index_set: the subset of indices that will be removed from the updated distance matrix
        @return: an updated distance matrix
        """
        if len(D) < self.large_matrix_size:
            return BuildTreeTopology.update_nj(D, index_set)
        else:
            return BuildTreeTopology.update_using_laplacian(D, index_set)


class Builder:
    """
    Instances store results as the program runs.
    """

    def __init__(self):
        self.scatter_points = []

    def get_R_lines(self, use_nj, use_modified_nj, use_all_spectral, use_one_spectral):
        """
        @return: a sequence of lines defining the header and data of an R table
        """
        # define the header line
        header_words = []
        header_words.append('atteson')
        if use_nj:
            header_words.append('nj')
        if use_modified_nj:
            header_words.append('nj.modified')
        if use_all_spectral:
            header_words.append('spectral.all')
        if use_one_spectral:
            header_words.append('spectral.one')
        header_line = '\t'.join(header_words)
        # define the data lines
        data_lines = []
        for i, p in enumerate(self.scatter_points):
            data_words = []
            data_words.append('%d' % (i+1))
            data_words.append('T' if p.atteson else 'F')
            if use_nj:
                data_words.append('%f' % p.nj_error)
            if use_modified_nj:
                data_words.append('%f' % p.modified_nj_error)
            if use_all_spectral:
                data_words.append('%f' % p.all_spectral_error)
            if use_one_spectral:
                data_words.append('%f' % p.one_spectral_error)
            data_line = '\t'.join(data_words)
            data_lines.append(data_line)
        # return lines of the R file
        return [header_line] + data_lines

    def evaluate(self, true_splits, D_estimated, atteson, use_nj, use_modified_nj, use_all_spectral, use_one_spectral):
        """
        @param true_splits: the set of all full splits implied by the true tree
        @param D_estimated: the estimated distance matrix
        @param atteson: True iff the distance matrix is Atteson
        """
        # initialize the errors
        nj_error = None
        modified_nj_error = None
        all_spectral_error = None
        one_spectral_error = None
        if use_nj:
            nj_splits = BuildTreeTopology.get_splits(D_estimated, BuildTreeTopology.split_nj, BuildTreeTopology.update_nj)
            nj_error = Xtree.splits_to_rf_distance(nj_splits, true_splits)
        if use_modified_nj:
            modified_nj_splits = BuildTreeTopology.get_splits(D_estimated, BuildTreeTopology.split_nj, BuildTreeTopology.update_using_laplacian)
            modified_nj_error = Xtree.splits_to_rf_distance(modified_nj_splits, true_splits)
        if use_all_spectral:
            splitter = BuildTreeTopology.split_using_eigenvector_with_nj_fallback
            updater = BuildTreeTopology.update_using_laplacian
            all_spectral_splits = BuildTreeTopology.get_splits(D_estimated, splitter, updater)
            all_spectral_error = Xtree.splits_to_rf_distance(all_spectral_splits, true_splits)
        if use_one_spectral:
            splitter = SplitFunctor(len(D_estimated))
            updater = UpdateFunctor(len(D_estimated))
            one_spectral_splits = BuildTreeTopology.get_splits(D_estimated, splitter, updater)
            one_spectral_error = Xtree.splits_to_rf_distance(one_spectral_splits, true_splits)
        # add the data point
        self.scatter_points.append(ScatterPoint(atteson, nj_error, modified_nj_error, all_spectral_error, one_spectral_error))


def process(ntaxa, length, nseconds, branch_length_sampler, use_nj, use_modified_nj, use_all_spectral, use_one_spectral):
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
                # determine whether or not the distance matrix is Atteson with respect to the tree
                atteson = BuildTreeTopology.is_atteson(tree, D)
                # record information about the splits
                builder.evaluate(true_splits, D, atteson, use_nj, use_modified_nj, use_all_spectral, use_one_spectral)
            except InfiniteDistanceError, e:
                n_infinite_errors += 1
            except ZeroDistanceError, e:
                n_zero_errors += 1
            except BuildTreeTopology.InvalidSpectralSplitException, e:
                n_failed_spectral_splits += 1
    except KeyboardInterrupt, e:
        pass
    # make the response
    out = StringIO()
    print >> out, '#', elapsed_time, 'seconds of run time'
    print >> out, '#', ntaxa, 'taxa per tree'
    print >> out, '#', branch_length_sampler
    print >> out, '#', length, 'nucleotides per sequence'
    print >> out, '#', n_zero_errors, 'samples were rejected because of an estimated distance of zero'
    print >> out, '#', n_infinite_errors, 'samples were rejected because of an estimated distance of infinity'
    print >> out, '#', n_failed_spectral_splits, 'distance matrices induced an intermediate distance matrix with no spectral split'
    for s in builder.get_R_lines(use_nj, use_modified_nj, use_all_spectral, use_one_spectral):
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
    response_text = process(ntaxa, length, nseconds, branch_length_sampler, fs.nj, fs.modified_nj, fs.all_spectral, fs.one_spectral)
    response_headers = [('Content-Type', 'text/plain')]
    if fs.attachment:
        output_filename = 'distances.table'
        response_headers.append(('Content-Disposition', "%s; filename=%s" % (fs.delivery, output_filename)))
    return response_headers, response_text

def main(options):
    assert 0 <= options.nseconds
    assert 4 <= options.ntaxa <= 20
    assert 10 <= options.length <= 10000
    branch_length_sampler = BranchLengthSampler.Pachter()
    print process(options.ntaxa, options.length, options.nseconds, branch_length_sampler, True, True, True, True)

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('--nseconds', dest='nseconds', type='int', default=0, help='seconds to run or 0 to run until ctrl-c')
    parser.add_option('--ntaxa', dest='ntaxa', type='int', default=20, help='number of taxa in each sampled tree topology')
    parser.add_option('--length', dest='length', type='int', default=200, help='nucleotide sequence length for distance matrix sampling')
    options, args = parser.parse_args()
    main(options)
