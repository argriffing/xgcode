"""Check the balancedness and quality of first spectral splits.

Check the balancedness and quality
of first spectral splits for a research paper.
Set a sequence length; for example 200 base pairs.
Sample some number of trees,
and from each tree sample an alignment and then a distance matrix.
For each distance matrix make an initial split,
and note the balancedness of the split and
whether or not it was compatible with the tree.
Return the results as an R table.
"""

from StringIO import StringIO
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
import Progress
import Form
import FormOut

class InfiniteDistanceError(Exception): pass
class ZeroDistanceError(Exception): pass

class TimeoutError(Exception): pass

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.Integer('ntaxa',
                'number of taxa for distance matrix sampling',
                20, low=4, high=20),
            Form.Integer('seqlen',
                'sequence length for distance matrix sampling',
                200, low=100, high=10000),
            Form.Integer('nsamples',
                'number of samples',
                10, low=1, high=1000),
            Form.RadioGroup('tree_sampling', 'branch length distribution', [
                Form.RadioItem('pachter_length',
                    str(BranchLengthSampler.Pachter()), True),
                Form.RadioItem('exponential_length',
                    str(BranchLengthSampler.Exponential())),
                Form.RadioItem('uniform_length_a',
                    str(BranchLengthSampler.UniformA())),
                Form.RadioItem('uniform_length_b',
                    str(BranchLengthSampler.UniformB()))]),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.RTable()

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


def get_response_content(fs):
    # allow only two seconds for web access
    nseconds = 2
    # read the options
    ntaxa = fs.ntaxa
    seqlen = fs.seqlen
    nsamples = fs.nsamples
    # define the branch length sampler
    if fs.pachter_length:
        branch_length_sampler = BranchLengthSampler.Pachter()
    elif fs.exponential_length:
        branch_length_sampler = BranchLengthSampler.Exponential()
    elif fs.uniform_length_a:
        branch_length_sampler = BranchLengthSampler.UniformA()
    elif fs.uniform_length_b:
        branch_length_sampler = BranchLengthSampler.UniformB()
    # return the response
    response_text = process(ntaxa, nseconds, seqlen, nsamples,
            branch_length_sampler, False)
    return response_text + '\n'

def process(ntaxa, nseconds, seqlen, nsamples, branch_length_sampler, use_pbar):
    """
    @param ntaxa: the number of taxa per tree
    @param nseconds: stop after this many seconds
    @param seqlen: use this sequence length
    @param nsamples: stop after this many samples per sequence length
    @param branch_length_sampler: this function samples branch lengths independently
    @param use_pbar: True iff a progress bar should be used
    @return: a multi-line string of the contents of an R table
    """
    # initialize the global rejection counts
    nrejected_zero = 0
    nrejected_inf = 0
    nrejected_fail = 0
    naccepted = 0
    # Initialize the accumulation matrix.
    # The rows specify the size of the smaller side of the initial split.
    # The columns specify the compatibility status of the split.
    nsmall_sizes = (ntaxa / 2) + 1
    accum = np.zeros((nsmall_sizes, 2), dtype=np.int)
    # Repeatedly analyze samples.
    # We might have to stop early if we run out of time or if ctrl-c is pressed.
    # If we have to stop early, then show the results of the progress so far.
    termination_reason = 'no reason for termination was given'
    start_time = time.time()
    pbar = Progress.Bar(nsamples) if use_pbar else None
    try:
        for sample_index in range(nsamples):
            # keep trying to get an accepted sample
            while True:
                # check the time
                if nseconds and time.time() - start_time > nseconds:
                    raise TimeoutError()
                # first sample a tree and get its set of informative splits
                tree = TreeSampler.sample_agglomerated_tree(ntaxa)
                true_splits = tree.get_nontrivial_splits()
                # sample the branch lengths
                for branch in tree.get_branches():
                    branch.length = branch_length_sampler()
                # Attempt to sample a distance matrix.
                # If the sample was rejected then note the reason and go back to the drawing board.
                try:
                    D = sample_distance_matrix(tree, seqlen)
                except InfiniteDistanceError as e
                    nrejected_inf += 1
                    continue
                except ZeroDistanceError as e
                    nrejected_zero += 1
                    continue
                # Attempt to estimate the primary split of the tree from the distance matrix.
                # If there was a technical failure then note it and go back to the drawing board.
                # Otherwise note the compatibility and balance of the split.
                try:
                    eigensplit = BuildTreeTopology.split_using_eigenvector(D)
                    small_size = min(len(side) for side in eigensplit)
                    if eigensplit in true_splits:
                        compatibility = 1
                    else:
                        compatibility = 0
                except BuildTreeTopology.DegenerateSplitException, e:
                    small_size = 0
                    compatibility = 1
                except BuildTreeTopology.InvalidSpectralSplitException, e:
                    nrejected_fail += 1
                    continue
                # reaching this point means the sample was accepted
                naccepted += 1
                break
            # accumulate the results of the accepted sample
            accum[small_size, compatibility] += 1
            if pbar:
                pbar.update(naccepted)
        else:
            termination_reason = 'the requested number of samples per sequence length was attained'
    except KeyboardInterrupt, e:
        termination_reason = 'keyboard interrupt'
    except TimeoutError as e
        termination_reason = 'time limit expired'
    if pbar:
        pbar.finish()
    # define the time limit string
    if nseconds:
        time_limit_string = 'the simulation was limited to %d seconds' % nseconds
    else:
        time_limit_string = 'no time limit was imposed'
    # create the results in convenient R table form
    out = StringIO()
    print >> out, '#', "R usage: mytable <- read.table('this.filename', header=T)"
    print >> out, '#', time.time() - start_time, 'elapsed seconds'
    print >> out, '#', time_limit_string
    print >> out, '#', 'the simulation was limited to', nsamples, 'samples per sequence length'
    print >> out, '#', 'reason for termination:', termination_reason
    print >> out, '#', seqlen, 'base pairs per sampled sequence'
    print >> out, '#', ntaxa, 'taxa per tree'
    print >> out, '#', branch_length_sampler
    print >> out, '#', naccepted, 'samples were accepted'
    print >> out, '#', nrejected_zero, 'samples were rejected with a distance estimate of zero'
    print >> out, '#', nrejected_inf, 'samples were rejected with a distance estimate of infinity'
    print >> out, '#', nrejected_fail, 'samples were rejected because of a technical problem with the split'
    print >> out, '\t'.join(('split.balance', 'nsplits.invalid', 'nsplits.valid'))
    for i, row in enumerate(accum):
        table_row = [str(i)] + [str(x) for x in row]
        print >> out, '\t'.join(table_row)
    # return the results
    return out.getvalue().strip()

def main(options):
    # validate the options
    assert 0 <= options.nseconds
    assert 4 <= options.ntaxa <= 30
    assert 1 <= options.seqlen
    assert 1 <= options.nsamples
    branch_length_sampler = BranchLengthSampler.UniformB()
    #branch_length_sampler = BranchLengthSampler.Pachter()
    use_pbar = True
    try:
        print process(options.ntaxa, options.nseconds, options.seqlen, options.nsamples, branch_length_sampler, use_pbar)
    except HandlingError as e
        print 'error:', e

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('--ntaxa', dest='ntaxa', type='int', default=20, help='number of taxa in each sampled tree topology')
    parser.add_option('--nseconds', dest='nseconds', type='int', default=0, help='seconds to run or 0 to run until ctrl-c')
    parser.add_option('--seqlen', dest='seqlen', type='int', default=200, help='sequence length')
    parser.add_option('--nsamples', dest='nsamples', type='int', default=10, help='number of samples')
    options, args = parser.parse_args()
    main(options)

