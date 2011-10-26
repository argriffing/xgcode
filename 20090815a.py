"""Simulate positions in aligned reads. [UNFINISHED]

The transition matrix is among the ordered hidden states
(homozygous (1), heterozygous (2), overcovered (x)).
Each row of the output corresponds to overlapping reads at a position.
For each row, the number of A, C, G, and T reads is given.
Optionally, the hidden state may be shown as well.
"""

from StringIO import StringIO
import time
import math
import random
import optparse

import numpy as np

from SnippetUtil import HandlingError
import MatrixUtil
import Form
import FormOut
import Progress
import FelTree

"""
def get_form():
    form_objects = [
            Form.Integer('ntaxa', 'number of taxa per tree', 20, low=4, high=20),
            Form.Integer('nsamples', 'number of trees to sample', 100, low=1, high=1000),
            Form.RadioGroup('tree_sampling', 'branch length distribution', [
                Form.RadioItem('pachter_length', str(BranchLengthSampler.Pachter()), True),
                Form.RadioItem('exponential_length', str(BranchLengthSampler.Exponential())),
                Form.RadioItem('uniform_length_a', str(BranchLengthSampler.UniformA())),
                Form.RadioItem('uniform_length_b', str(BranchLengthSampler.UniformB()))])]
    return form_objects

    parser.add_option('--ntaxa', dest='ntaxa', type='int', default=20, help='number of taxa in each sampled tree topology')
    parser.add_option('--nseconds', dest='nseconds', type='int', default=0, help='seconds to run or 0 to run until ctrl-c')
    parser.add_option('--nsamples', dest='nsamples', type='int', default=100, help='number of samples')

    assert 0 <= options.nseconds
    assert 4 <= options.ntaxa <= 20
    assert 1 <= options.nsamples
    branch_length_sampler = BranchLengthSampler.UniformB()
    use_pbar = True
    print process(options.ntaxa, options.nseconds, options.nsamples, branch_length_sampler, use_pbar)
"""


def get_form():
    """
    @return: the body of a form
    """
    # FIXME
    T = []
    # define the default transition matrix
    # define the objects
    form_objects = [
            Form.Integer('good_coverage',
                'expected coverage of homozygous and heterozygous positions',
                10, low=1),
            Form.Integer('bad_coverage',
                'expected coverage of overcovered positions',
                30, low=1),
            Form.Float('randomization_rate',
                'read randomization rate',
                0.1, low_inclusive=0, high_inclusive=1),
            Form.Integer('npositions',
                'number of positions to sample',
                100, low=1, high=10000),
            Form.Matrix('transition_matrix',
                'transition matrix',
                T, MatrixUtil.assert_transition_matrix),
            Form.CheckGroup('options', 'output options', [
                Form.CheckItem('show_hidden_states',
                    'show hidden states', True)]),
            Form.Integer('seed',
                'prng seed (or a negative for no specified seed)', -1)]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):
    # allow only two seconds for web access
    nseconds = 2
    # read the options
    ntaxa = fs.ntaxa
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
    # get the response
    response_text = process(ntaxa, nseconds, nsamples,
            branch_length_sampler, False)
    return response_text + '\n'

def sample_perturbation_matrix(n, frobnorm):
    """
    The sampled matrix is symmetric and has zeros on the diagonal.
    @param n: the number of rows (and columns) in the sampled matrix
    @param frobnorm: the Frobenius norm of the sampled matrix
    @return: an nxn numpy array with a given Frobenius norm
    """
    # initialize the matrix to all zeros
    E = np.zeros((n,n))
    # define the symmetric off-diagonal perturbations
    for i in range(1, n):
        for j in range(i+1, n):
            value = random.gauss(0, 1)
            E[i, j] = value
            E[j, i] = value
    # force the matrix to have a given Frobenius norm
    return frobnorm * (E / np.linalg.norm(E))

def get_sorted_eigensystem(M):
    w, vT = np.linalg.eigh(M)
    w_v_pairs = zip(w, vT.T.tolist())
    return list(sorted(w_v_pairs))

def get_stability(D):
    """
    The stability is defined as a bound on norms of perturbation matrices.
    If D is perturbed by a matrix whose Frobenius norm is less than the stability,
    then the spectral split remains unchanged.
    @param D: a distance matrix
    @return: the stability of the distance matrix
    """
    HDH = MatrixUtil.double_centered(D)
    # get the eigendecomposition
    w_v_pairs = get_sorted_eigensystem(-HDH)
    # compute the eigengap
    w = [w for w, v in w_v_pairs]
    lambda_1 = w[-1]
    lambda_2 = w[-2]
    eigengap = lambda_1 - lambda_2
    delta = eigengap
    # compute an eigenvector stability statistic
    v = [v for w, v in w_v_pairs]
    dominant_eigenvector = v[-1]
    alpha = min(abs(x) for x in dominant_eigenvector)
    # compute the stability as a function of alpha and delta
    eigenvalue_control = delta / (2*math.sqrt(2))
    eigenvector_control = alpha * delta / (4 + math.sqrt(2)*alpha) 
    stability = min(eigenvalue_control, eigenvector_control)
    return stability

def get_split(D):
    """
    @param D: an exact or perturbed distance matrix
    @return: a set of frozensets of indices
    """
    HDH = MatrixUtil.double_centered(D)
    # get the dominant eigenvector
    w_v_pairs = get_sorted_eigensystem(-HDH)
    v = [v for w, v in w_v_pairs]
    ev_dom = v[-1]
    neg_set = frozenset(i for i, x in enumerate(ev_dom) if x < 0)
    nonneg_set = frozenset(i for i, x in enumerate(ev_dom) if x >= 0)
    return set([neg_set, nonneg_set])

def process(ntaxa, nseconds, nsamples, branch_length_sampler, use_pbar):
    """
    @param ntaxa: the number of taxa per tree
    @param nseconds: stop after this many seconds
    @param nsamples: stop after this many samples
    @param branch_length_sampler: this function samples branch lengths independently
    @param use_pbar: True iff a progress bar should be used
    @return: a multi-line string of the contents of an R table
    """
    a_successes = 0
    a_failures = 0
    b_successes = 0
    b_failures = 0
    # Repeatedly analyze samples.
    # We might have to stop early if we run out of time or if ctrl-c is pressed.
    # If we have to stop early, then show the results of the progress so far.
    termination_reason = 'no reason for termination was given'
    start_time = time.time()
    pbar = Progress.Bar(nsamples) if use_pbar else None
    try:
        for sample_index in range(nsamples):
            # check the time
            if nseconds and time.time() - start_time > nseconds:
                raise TimeoutError()
            # sample a tree
            tree = TreeSampler.sample_agglomerated_tree(ntaxa)
            for branch in tree.get_branches():
                branch.length = branch_length_sampler()
            D = np.array(tree.get_distance_matrix())
            # get the split defined by the tree
            original_split = get_split(D)
            # get the stability of the split
            stability = get_stability(D)
            # sample a perturbation matrix that should not change the split
            E = sample_perturbation_matrix(ntaxa, stability/2)
            # evaluate the split induced by the unerperturbed perturbed distance matrix
            perturbed_split = get_split(D + E)
            if original_split == perturbed_split:
                a_successes += 1
            else:
                a_failures += 1
            # evaluage the split induced by the overperturbed distance matrix
            perturbed_split = get_split(D + E*200)
            if original_split == perturbed_split:
                b_successes += 1
            else:
                b_failures += 1
            # update the progress bar
            if pbar:
                pbar.update(sample_index + 1)
        else:
            termination_reason = 'the requested number of samples was attained'
    except KeyboardInterrupt, e:
        termination_reason = 'keyboard interrupt'
    except TimeoutError, e:
        termination_reason = 'time limit expired'
    if pbar:
        pbar.finish()
    # define the time limit string
    if nseconds:
        time_limit_string = 'the simulation was limited to %d seconds' % nseconds
    else:
        time_limit_string = 'no time limit was imposed'
    # create the results
    out = StringIO()
    print >> out, 'debug information:'
    print >> out, time.time() - start_time, 'elapsed seconds'
    print >> out, time_limit_string
    print >> out, 'the simulation was limited to', nsamples, 'samples'
    print >> out, 'reason for termination:', termination_reason
    print >> out, ntaxa, 'taxa per tree'
    print >> out, branch_length_sampler
    print >> out
    print >> out, 'results:'
    print >> out, a_successes, '\tsuccesses when success is guaranteed'
    print >> out, a_failures, '\tfailures when success is guaranteed'
    print >> out, b_successes, '\tsuccesses when success is not guaranteed'
    print >> out, b_failures, '\tfailures when success is not guaranteed'
    # return the results
    return out.getvalue().strip()

def main(options):
    # validate the options
    assert 0 <= options.nseconds
    assert 4 <= options.ntaxa <= 20
    assert 1 <= options.nsamples
    branch_length_sampler = BranchLengthSampler.UniformB()
    use_pbar = True
    print process(options.ntaxa, options.nseconds, options.nsamples, branch_length_sampler, use_pbar)

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('--ntaxa', dest='ntaxa', type='int', default=20, help='number of taxa in each sampled tree topology')
    parser.add_option('--nseconds', dest='nseconds', type='int', default=0, help='seconds to run or 0 to run until ctrl-c')
    parser.add_option('--nsamples', dest='nsamples', type='int', default=100, help='number of samples')
    options, args = parser.parse_args()
    main(options)

