"""Compare a spectral tree topology estimation method to neighbor joining.

Several phylogenetic tree topologies
and corresponding approximate distance matrices are sampled,
and for each distance matrix a tree topology
is estimated using each of two methods.
Phylogenetic trees are sampled by random agglomeration of some number of taxa.
A distance matrix is sampled given a phylogenetic tree by first sampling
a nucleotide alignment according to a Jukes-Cantor model,
and then estimating the pairwise distances
by the maximum likelihood Jukes-Cantor distance.
A sample is rejected when a pairwise distance is zero or infinity.
"""

from StringIO import StringIO
import time
import math
import random

import numpy as np

from SnippetUtil import HandlingError
import MatrixUtil
import JC69
import TreeSampler
import BuildTreeTopology
import BranchLengthSampler
import Form
import FormOut

class InfiniteDistanceError(Exception): pass
class ZeroDistanceError(Exception): pass

#FIXME use the distance matrix sampling library

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.Integer('length', 'sequence length',
                1000, low=10, high=20000),
            Form.Integer('ntaxa', 'number of taxa',
                20, low=4, high=20),
            Form.RadioGroup('tree_sampling', 'branch length distribution', [
                Form.RadioItem('pachter_length',
                    str(BranchLengthSampler.Pachter()), True),
                Form.RadioItem('exponential_length',
                    str(BranchLengthSampler.Exponential())),
                Form.RadioItem('uniform_length_a',
                    str(BranchLengthSampler.UniformA())),
                Form.RadioItem('uniform_length_b',
                    str(BranchLengthSampler.UniformB()))]),
            Form.RadioGroup('first_method', 'first estimation method', [
                Form.RadioItem('first_nj', 'neighbor joining', True),
                Form.RadioItem('first_modnj', 'modified neighbor joining'),
                Form.RadioItem('first_specnj', 'spectral reconstruction')]),
            Form.RadioGroup('second_method', 'second estimation method', [
                Form.RadioItem('second_nj', 'neighbor joining'),
                Form.RadioItem('second_modnj', 'modified neighbor joining'),
                Form.RadioItem('second_specnj', 'spectral', True)])]
    return form_objects

def get_form_out():
    return FormOut.Report()

def sample_distance_matrix(xtree_root, sequence_length):
    sequences = JC69.sample_xtree_sequences(xtree_root, sequence_length)
    nsequences = len(sequences)
    pairwise_mismatch_count = np.zeros((nsequences, nsequences))
    for i, sa in enumerate(sequences):
        for j, sb in enumerate(sequences):
            if i < j:
                nmismatches = sum(1 for a, b in zip(sa, sb) if a != b)
                if not nmismatches:
                    raise ZeroDistanceError()
                if nmismatches * 4 >= sequence_length * 3:
                    raise InfiniteDistanceError()
                pairwise_mismatch_count[i][j] = nmismatches
    D = np.zeros_like(pairwise_mismatch_count)
    for i in range(nsequences):
        for j in range(nsequences):
            if i < j:
                d_raw = pairwise_mismatch_count[i][j] / float(sequence_length)
                b = 0.75
                d_mle = -b*math.log(1 - d_raw/b)
                D[i][j] = d_mle
                D[j][i] = d_mle
    return D


class Builder:

    def __init__(self, split_function, update_function, name):
        self.split_function = split_function
        self.update_function = update_function
        self.name = name

    def evaluate(self, true_splits, D_estimated):
        """
        @param true_splits: a set of full splits that defines the true tree topology
        @param D_estimated: an estimated distance matrix conformant to the split labels
        @return: 1 if success, 0 if failure
        """
        estimated_splits = BuildTreeTopology.get_splits(D_estimated, self.split_function, self.update_function)
        if estimated_splits == true_splits:
            return 1
        else:
            return 0


def pachter_branch_length_sampler():
    """
    This length is 0.1 as in the sampling procedure in the paper "Why neighbor joining works".
    """
    return 0.1

def uniform_branch_length_sampler():
    """
    Each length is chosen according to a uniform (0, 1) distribution.
    """
    return random.random()

def process(ntaxa, length, nseconds, builders, branch_length_sampler):
    """
    @param ntaxa: the number of taxa in the sampled trees
    @param length: the length of sequences used to sample the distance matrix
    @param nseconds: allow this many seconds to run
    @param builders: tree builder objects
    @param branch_length_sampler: returns a tree drawn from some distribution
    @return: a multi-line string that summarizes the results
    """
    start_time = time.time()
    # track the number of samples that failed for various reasons
    n_zero_errors = 0
    n_infinite_errors = 0
    n_failed_spectral_splits = 0
    # define the number of attempts that fall into each of the four categories
    non_atteson_results = [[0, 0], [0, 0]]
    atteson_results = [[0, 0], [0, 0]]
    #pachter_results = [[0, 0], [0, 0]]
    # evaluate the quality of reconstructions from a bunch of different samples
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
                a, b = [builder.evaluate(true_splits, D) for builder in builders]
                if BuildTreeTopology.is_atteson(tree, D):
                    atteson_results[a][b] += 1
                #elif BuildTreeTopology.is_quartet_additive(tree, D) and BuildTreeTopology.is_quartet_consistent(tree, D):
                    #pachter_results[a][b] += 1
                else:
                    non_atteson_results[a][b] += 1
            except InfiniteDistanceError as e:
                n_infinite_errors += 1
            except ZeroDistanceError as e:
                n_zero_errors += 1
            except BuildTreeTopology.InvalidSpectralSplitException, e:
                n_failed_spectral_splits += 1
    except KeyboardInterrupt, e:
        pass
    # make the response
    result_to_string = ['failed', 'succeeded']
    out = StringIO()
    print >> out, elapsed_time, 'seconds spent on this group'
    print >> out, n_zero_errors, 'samples were rejected because of an estimated distance of zero'
    print >> out, n_infinite_errors, 'samples were rejected because of an estimated distance of infinity'
    print >> out, n_failed_spectral_splits, 'distance matrices induced an intermediate distance matrix with no spectral split'
    print >> out
    print >> out, 'results for distance matrices that satisfy the atteson condition:'
    for first_result in (0, 1):
        sa = result_to_string[first_result] + ' for ' + builders[0].name
        for second_result in (0, 1):
            sb = result_to_string[second_result] + ' for ' + builders[1].name
            nresults = atteson_results[first_result][second_result]
            print >> out, nresults, 'reconstructions', sa, 'and', sb
    print >> out
    """
    print >> out, 'results for remaining distance matrices that are quartet additive and consistent:'
    for first_result in (0, 1):
        sa = result_to_string[first_result] + ' for ' + builders[0].name
        for second_result in (0, 1):
            sb = result_to_string[second_result] + ' for ' + builders[1].name
            nresults = pachter_results[first_result][second_result]
            print >> out, nresults, 'reconstructions', sa, 'and', sb
    print >> out
    """
    print >> out, 'results for the remaining distance matrices:'
    for first_result in (0, 1):
        sa = result_to_string[first_result] + ' for ' + builders[0].name
        for second_result in (0, 1):
            sb = result_to_string[second_result] + ' for ' + builders[1].name
            nresults = non_atteson_results[first_result][second_result]
            print >> out, nresults, 'reconstructions', sa, 'and', sb
    return out.getvalue().strip()

def get_response_content(fs):
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
    # define the first builder
    if fs.first_nj:
        splitter = BuildTreeTopology.split_nj
        updater = BuildTreeTopology.update_nj
        first_builder = Builder(splitter, updater, 'nj')
    elif fs.first_modnj:
        splitter = BuildTreeTopology.split_nj
        updater = BuildTreeTopology.update_using_laplacian
        first_builder = Builder(splitter, updater, 'modnj')
    elif fs.first_specnj:
        splitter = BuildTreeTopology.split_using_eigenvector_with_nj_fallback
        updater = BuildTreeTopology.update_using_laplacian
        first_builder = Builder(splitter, updater, 'specnj')
    # define the second builder
    if fs.second_nj:
        splitter = BuildTreeTopology.split_nj
        updater = BuildTreeTopology.update_nj
        second_builder = Builder(splitter, updater, 'nj')
    elif fs.second_modnj:
        splitter = BuildTreeTopology.split_nj
        updater = BuildTreeTopology.update_using_laplacian
        second_builder = Builder(splitter, updater, 'modnj')
    elif fs.second_specnj:
        splitter = BuildTreeTopology.split_using_eigenvector_with_nj_fallback
        updater = BuildTreeTopology.update_using_laplacian
        second_builder = Builder(splitter, updater, 'specnj')
    builders = [first_builder, second_builder]
    response_text = process(
            ntaxa, length, nseconds, builders, branch_length_sampler)
    return response_text

def main():
    ntaxa = 20
    length = 200
    nseconds = 5*60
    print ntaxa, 'taxa per tree'
    print length, 'nucleotides per sequence'
    print
    branch_length_samplers = (BranchLengthSampler.Pachter(), BranchLengthSampler.UniformB())
    for i, branch_length_sampler in enumerate(branch_length_samplers):
        print 'group', i+1
        print branch_length_sampler
        print
        builders = [
            Builder(BuildTreeTopology.split_nj, BuildTreeTopology.update_nj, 'nj'),
            Builder(BuildTreeTopology.split_using_eigenvector_with_nj_fallback, BuildTreeTopology.update_using_laplacian, 'specnj')]
        print process(ntaxa, length, nseconds, builders, branch_length_sampler)
        print

if __name__ == '__main__':
    main()
