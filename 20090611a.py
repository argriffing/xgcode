"""Find an Atteson distance matrix and tree for which the first eigensplit fails.

The first eigensplit is the split defined
by the signs of the loadings of the first principal eigenvector of
the doubly centered distance matrix.
This split is considered to have failed if
there is a branch in the tree such that
each of the two subtrees defined by the branch
have at least one positively valuated tip
and one negatively valuated tip.
"""

import StringIO
import time
import random
import optparse

import numpy

from SnippetUtil import HandlingError
import MatrixUtil
import BuildTreeTopology
import BranchLengthSampler
import TreeSampler
import Xtree
import Form

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.Integer('ntaxa', 'number of taxa', 20, low=4, high=20),
            Form.RadioGroup('tree_sampling', 'branch length distribution', [
                Form.RadioItem('pachter_length', str(BranchLengthSampler.Pachter()), True),
                Form.RadioItem('exponential_length', str(BranchLengthSampler.Exponential())),
                Form.RadioItem('uniform_length_a', str(BranchLengthSampler.UniformA())),
                Form.RadioItem('uniform_length_b', str(BranchLengthSampler.UniformB()))])]
    return form_objects

def sample_atteson_distance_matrix(xtree_root):
    """
    @param xtree_root: the root of a weighted phylogenetic xtree
    """
    # get the minimum branch length
    shortest_branch_length = min(branch.length for branch in xtree_root.get_branches())
    # get the true distance matrix
    D = numpy.array(xtree_root.get_distance_matrix())
    n = len(D)
    # perturb entries of the distance matrix in a way that satisfies the Atteson condition
    for i in range(n):
        for j in range(n):
            if i < j:
                exact_distance = D[i][j]
                noise = random.uniform(-0.5, 0.5) * shortest_branch_length
                perturbed_distance = exact_distance + noise
                D[i][j] = perturbed_distance
                D[j][i] = perturbed_distance
    return D

def process(ntaxa, nseconds, branch_length_sampler):
    """
    @param ntaxa: the number of taxa in the sampled trees
    @param nseconds: allow this many seconds to run or None to run forever
    @param branch_length_sampler: a functor that returns a branch length and has a string cast
    @return: a multi-line string that summarizes the results
    """
    start_time = time.time()
    # initialize some state that will be tracked over the entire run
    degenerate_count = 0
    invalid_split_count = 0
    valid_split_count = 0
    spectral_error_count = 0
    atteson_error_count = 0
    counterexample_D = None
    counterexample_tree = None
    # do a bunch of reconstructions from sampled distance matrices
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
            # sample the atteson distance matrix
            D = sample_atteson_distance_matrix(tree)
            # assert that the atteson condition is true
            if not BuildTreeTopology.is_atteson(tree, D):
                atteson_error_count += 1
            else:
                try:
                    # see if the eigensplit is in the set of true splits
                    eigensplit = BuildTreeTopology.split_using_eigenvector(D)
                    if eigensplit in true_splits:
                        valid_split_count += 1
                    else:
                        invalid_split_count += 1
                        counterexample_D = D
                        counterexample_tree = tree
                        break
                except BuildTreeTopology.DegenerateSplitException, e:
                    degenerate_count += 1
                except BuildTreeTopology.InvalidSpectralSplitException, e:
                    spectral_error_count += 1
    except KeyboardInterrupt, e:
        pass
    # make the response
    out = StringIO.StringIO()
    if invalid_split_count:
        print >> out, 'A counterexample was found!'
        print >> out, 'D:'
        print >> out, MatrixUtil.m_to_string(counterexample_D)
        print >> out, 'tree:'
        print >> out, counterexample_tree.get_newick_string()
    else:
        print >> out, 'No counterexample was found.'
    print >> out, elapsed_time, 'seconds of run time'
    print >> out, ntaxa, 'taxa per tree'
    print >> out, branch_length_sampler
    print >> out, spectral_error_count, 'spectral errors (should never happen)'
    print >> out, atteson_error_count , 'atteson errors (should never happen)'
    print >> out, degenerate_count, 'degenerate splits found'
    print >> out, invalid_split_count, 'invalid splits found (should be zero or one)'
    print >> out, valid_split_count, 'valid splits found'
    return out.getvalue().strip()

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    nseconds = 2
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
    response_text = process(ntaxa, nseconds, branch_length_sampler)
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, response_text

def main(options):
    assert 0 <= options.nseconds
    assert 4 <= options.ntaxa <= 20
    branch_length_sampler = BranchLengthSampler.Pachter()
    print process(options.ntaxa, options.nseconds, branch_length_sampler)

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('--nseconds', dest='nseconds', type='int', default=0, help='seconds to run or 0 to run until ctrl-c')
    parser.add_option('--ntaxa', dest='ntaxa', type='int', default=20, help='number of taxa in each sampled tree topology')
    options, args = parser.parse_args()
    main(options)
