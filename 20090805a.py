"""Do a tree reconstruction simulation for the research paper.

For each of N sequence lengths, sample K trees;
from each tree sample a distance matrix;
for each distance matrix, reconstruct the tree using different methods.
"""

import StringIO
import time
import math
import random
import optparse

import numpy as np

from SnippetUtil import HandlingError
import MatrixUtil
import BranchLengthSampler
import TreeSampler
import BuildTreeTopology
import Xtree
import JC69
import Euclid
import Form
import Progress

# define the headers of the R table
g_headers = [
        'sequence.length',
        'nsamples.accepted',
        'nsamples.accepted.atteson',
        'nsamples.rejected.zero',
        'nsamples.rejected.inf',
        'nsamples.rejected.fail',
        'nsuccesses.both',
        'nsuccesses.neither',
        'nsuccesses.nj.only',
        'nsuccesses.topdown.only',
        'first.split.informative',
        'first.split.uninformative',
        'first.split.invalid']

class InfiniteDistanceError(Exception): pass
class ZeroDistanceError(Exception): pass

class TimeoutError(Exception): pass

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.Integer('ntaxa', 'number of taxa for distance matrix sampling', 20, low=4, high=20),
            Form.Integer('nlengths', 'number of sequence lengths to consider (must be odd)', 9, low=3, high=99),
            Form.Integer('nsamples', 'number of samples per sequence length', 5, low=1, high=99),
            Form.RadioGroup('tree_sampling', 'branch length distribution', [
                Form.RadioItem('pachter_length', str(BranchLengthSampler.Pachter()), True),
                Form.RadioItem('exponential_length', str(BranchLengthSampler.Exponential())),
                Form.RadioItem('uniform_length_a', str(BranchLengthSampler.UniformA())),
                Form.RadioItem('uniform_length_b', str(BranchLengthSampler.UniformB()))]),
            Form.RadioGroup('distance_options', 'recursive matrix construction method', [
                Form.RadioItem('pruning_like', 'go through the Laplacian, like Felsenstein pruning', True),
                Form.RadioItem('nj_like', 'directly use distances, like neighbor joining')]),
            Form.RadioGroup('delivery', 'delivery', [
                Form.RadioItem('inline', 'view as text', True),
                Form.RadioItem('attachment', 'download as an R table')])]
    return form_objects

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


def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # allow only two seconds for web access
    nseconds = 2
    # read the options
    ntaxa = fs.ntaxa
    nlengths = fs.nlengths
    nsamples = fs.nsamples
    nj_like = fs.nj_like
    # do extra validation
    if nlengths % 2 == 0:
        raise HandlingError('the number of sequence lengths must be odd')
    # define the branch length sampler
    if fs.pachter_length:
        branch_length_sampler = BranchLengthSampler.Pachter()
    elif fs.exponential_length:
        branch_length_sampler = BranchLengthSampler.Exponential()
    elif fs.uniform_length_a:
        branch_length_sampler = BranchLengthSampler.UniformA()
    elif fs.uniform_length_b:
        branch_length_sampler = BranchLengthSampler.UniformB()
    # get the response
    response_text = process(ntaxa, nseconds, nlengths, nsamples, nj_like, branch_length_sampler, False)
    response_headers = [('Content-Type', 'text/plain')]
    if fs.attachment:
        response_headers.append(('Content-Disposition', 'attachment; filename=distances.table'))
    return response_headers, response_text

def gen_sequence_lengths_helper(n, low, high):
    """
    Generate integer lengths that are nearly evenly spaced on a log scale.
    @param n: the number of sequence lengths to generate
    @param low: the smallest sequence length
    @param high: the highest sequence length
    """
    yield low
    incr = (float(high) / float(low)) ** (1.0 / (n-1))
    for i in range(1, n-1):
        yield int(low*(incr ** i))
    yield high

def get_sequence_lengths(nlengths):
    """
    Each length is a positive integer between 100 and 10000.
    The lengths 100, 1000, and 10000 are always included.
    The ratios of consecutive lengths are roughly constant.
    @param nlengths: the number of sequence lengths
    @return: a list of sequence lengths
    """
    # the number of lengths must be odd and at least 3
    assert nlengths % 2
    assert 3 <= nlengths
    first_lengths = list(gen_sequence_lengths_helper((nlengths+1)/2, 100, 1000))
    second_lengths = list(gen_sequence_lengths_helper((nlengths+1)/2, 1000, 10000))
    lengths = first_lengths + second_lengths[1:]
    return lengths

def incr_attribute(attribute_array, attribute):
    """
    Increment an element in an attribute array.
    The motivation of this approach is the convenience of numpy array addition.
    @param attribute_array: a numpy array conformant to the global header list
    @param attribute: an element of the global header list
    @return: the attribute array
    """
    header_to_index = dict((header, i) for i, header in enumerate(g_headers))
    index = header_to_index[attribute]
    attribute_array[index] += 1
    return attribute_array

def get_attribute(attribute_array, attribute):
    """
    Get an element in an attribute array.
    The motivation of this approach is the convenience of numpy array addition.
    @param attribute_array: an array conformant to the global header list
    @param attribute: an element of the global header list
    @return: an element of the attribute array
    """
    header_to_index = dict((header, i) for i, header in enumerate(g_headers))
    index = header_to_index[attribute]
    return attribute_array[index]

def set_attribute(attribute_array, attribute, value):
    """
    Set an element in an attribute array.
    The motivation of this approach is the convenience of numpy array addition.
    @param attribute_array: an array conformant to the global header list
    @param attribute: an element of the global header list
    @param value: the value to set
    """
    header_to_index = dict((header, i) for i, header in enumerate(g_headers))
    index = header_to_index[attribute]
    attribute_array[index] = value

def get_sample_results(sequence_length, ntaxa, nj_like, branch_length_sampler):
    """
    @param sequence_length: the length of each sequence in the sampled alignment
    @param ntaxa: the number of sequences in the sampled tree
    @param nj_like: True to create subsequent distance matrices using a generalized neighbor-joining-like approach
    @param branch_length_sampler: the length of each branch is independently sampled by this function
    @return: a numpy array conformant to the global header list
    """
    # initialize the array that will be returned
    attribute_array = np.zeros((len(g_headers),), dtype=np.int)
    # first sample a tree and get its set of informative splits
    tree = TreeSampler.sample_agglomerated_tree(ntaxa)
    true_splits = tree.get_nontrivial_splits()
    # sample the branch lengths
    for branch in tree.get_branches():
        branch.length = branch_length_sampler()
    # sample a distance matrix
    try:
        D = sample_distance_matrix(tree, sequence_length)
    except InfiniteDistanceError, e:
        return incr_attribute(attribute_array, 'nsamples.rejected.inf')
    except ZeroDistanceError, e:
        return incr_attribute(attribute_array, 'nsamples.rejected.zero')
    except BuildTreeTopology.InvalidSpectralSplitException, e:
        return incr_attribute(attribute_array, 'nsamples.rejected.fail')
    # see if the top down reconstruction was successful
    try:
        splitter = BuildTreeTopology.split_using_eigenvector_with_nj_fallback
        if nj_like:
            updater = BuildTreeTopology.update_generalized_nj
        else:
            updater = BuildTreeTopology.update_using_laplacian
        all_spectral_splits = BuildTreeTopology.get_splits(D, splitter, updater)
        top_down_success = (all_spectral_splits == true_splits)
    except BuildTreeTopology.InvalidSpectralSplitException, e:
        return incr_attribute(attribute_array, 'nsamples.rejected.fail')
    # at this point the sample is accepted
    incr_attribute(attribute_array, 'nsamples.accepted')
    # determine whether or not the distance matrix is Atteson with respect to the tree
    if BuildTreeTopology.is_atteson(tree, D):
        incr_attribute(attribute_array, 'nsamples.accepted.atteson')
    # see if the bottom up reconstruction was successful
    nj_splits = BuildTreeTopology.get_splits(D, BuildTreeTopology.split_nj, BuildTreeTopology.update_nj)
    nj_success = (nj_splits == true_splits)
    # note the joint results of the two reconstruction methods
    if top_down_success and nj_success:
        incr_attribute(attribute_array, 'nsuccesses.both')
    elif (not top_down_success) and (not nj_success):
        incr_attribute(attribute_array, 'nsuccesses.neither')
    elif top_down_success and (not nj_success):
        incr_attribute(attribute_array, 'nsuccesses.topdown.only')
    elif (not top_down_success) and nj_success:
        incr_attribute(attribute_array, 'nsuccesses.nj.only')
    # characterize the result of the first spectral split
    try:
        eigensplit = BuildTreeTopology.split_using_eigenvector(D)
        if eigensplit in true_splits:
            incr_attribute(attribute_array, 'first.split.informative')
        else:
            incr_attribute(attribute_array, 'first.split.invalid')
    except BuildTreeTopology.DegenerateSplitException, e:
        incr_attribute(attribute_array, 'first.split.uninformative')
    # return the attribute array
    return attribute_array

def process(ntaxa, nseconds, nlengths, nsamples, nj_like, branch_length_sampler, use_pbar):
    """
    @param ntaxa: the number of taxa per tree
    @param nseconds: stop after this many seconds
    @param nlengths: use this many different sequence lengths
    @param nsamples: stop after this many samples per sequence length
    @param nj_like: True to use a generalized neighbor-joining-like method of computing successive distance matrices
    @param branch_length_sampler: this function samples branch lengths independently
    @param use_pbar: True iff a progress bar should be used
    @return: a multi-line string of the contents of an R table
    """
    # define the sequence lengths
    lengths = get_sequence_lengths(nlengths)
    # initialize the accumulation matrix
    accum = np.zeros((nlengths, len(g_headers)), dtype=np.int)
    for i, sequence_length in enumerate(lengths):
        set_attribute(accum[i], 'sequence.length', sequence_length)
    # Repeatedly analyze samples from each sequence length.
    # We might have to stop early if we run out of time or if ctrl-c is pressed.
    # If we have to stop early, then show the results of the progress so far.
    termination_reason = 'no reason for termination was given'
    start_time = time.time()
    pbar = None
    if use_pbar:
        pbar = Progress.Bar(nsamples)
    try:
        for sample_index in range(nsamples):
            # reset the accumulation matrix for this iteration
            single_iteration_accum = np.zeros((nlengths, len(g_headers)))
            # accumulate attributes of sampling attempts for each sequence length
            for sequence_length_index, sequence_length in enumerate(lengths):
                # keep trying to get an accepted sample
                while True:
                    # check the time
                    if nseconds and time.time() - start_time > nseconds:
                        raise TimeoutError()
                    # get counts of attributes of a sample
                    sample_result = get_sample_results(sequence_length, ntaxa, nj_like, branch_length_sampler)
                    single_iteration_accum[sequence_length_index] += sample_result
                    # if the sample was accepted then we are done looking
                    if get_attribute(sample_result, 'nsamples.accepted'):
                        break
            # finish the iteration
            accum += single_iteration_accum
            if pbar:
                pbar.update(sample_index + 1)
        else:
            termination_reason = 'the requested number of samples per sequence length was attained'
    except KeyboardInterrupt, e:
        termination_reason = 'keyboard interrupt'
    except TimeoutError, e:
        termination_reason = 'time limit expired'
    if pbar:
        pbar.finish()
    # define the matrix successor creation method explanatory string
    if nj_like:
        matrix_successor_explanation = 'like neighbor joining'
    else:
        matrix_successor_explanation = 'like Felsenstein pruning'
    # define the time limit string
    if nseconds:
        time_limit_string = '%d seconds' % nseconds
    else:
        time_limit_string = '(no time limit)'
    # create the results in convenient R table form
    out = StringIO.StringIO()
    print >> out, '#', 'R usage: mytable <- read.table(\'this.filename\')'
    print >> out, '#', time.time() - start_time, 'elapsed seconds'
    print >> out, '#', 'the simulation was limited to', time_limit_string
    print >> out, '#', 'the simulation was limited to', nsamples, 'samples per sequence length'
    print >> out, '#', 'reason for termination:', termination_reason
    print >> out, '#', 'matrix successor creation method:', matrix_successor_explanation
    print >> out, '#', ntaxa, 'taxa per tree'
    print >> out, '#', branch_length_sampler
    print >> out, '\t'.join(g_headers)
    for i, row in enumerate(accum):
        print >> out, '\t'.join(str(x) for x in [i+1] + row.tolist())
    # return the results
    return out.getvalue().strip()

def main(options):
    # validate the options
    assert 0 <= options.nseconds
    assert 4 <= options.ntaxa <= 30
    assert options.nlengths % 2
    assert 3 <= options.nlengths
    assert 1 <= options.nsamples
    branch_length_sampler = BranchLengthSampler.UniformB()
    #branch_length_sampler = BranchLengthSampler.Pachter()
    use_pbar = True
    try:
        print process(options.ntaxa, options.nseconds, options.nlengths, options.nsamples, options.nj_like, branch_length_sampler, use_pbar)
    except HandlingError, e:
        print 'Error:', e

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('--ntaxa', dest='ntaxa', type='int', default=20, help='number of taxa in each sampled tree topology')
    parser.add_option('--nseconds', dest='nseconds', type='int', default=0, help='seconds to run or 0 to run until ctrl-c')
    parser.add_option('--nlengths', dest='nlengths', type='int', default=29, help='number of sequence lengths')
    parser.add_option('--nsamples', dest='nsamples', type='int', default=5, help='number of samples per sequence length')
    parser.add_option('--nj-like', action='store_true', dest='nj_like', default=False, help='use a generalized NJ-like way to create successor distance matrices')
    options, args = parser.parse_args()
    main(options)

