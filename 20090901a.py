"""Compare the Atteson and spectral bounds for some trees.
"""

from StringIO import StringIO
import time
import random
import optparse
import math

import numpy as np

from SnippetUtil import HandlingError
import MatrixUtil
import BranchLengthSampler
import TreeSampler
import Xtree
import Form
import FormOut

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.Integer('ntaxa', 'number of taxa', 20, low=4, high=20),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.RTable()


class BranchLengthComboPack:

    def __init__(self):
        self.components = [
            BranchLengthSampler.Pachter(),
            BranchLengthSampler.Exponential(),
            BranchLengthSampler.UniformA(),
            BranchLengthSampler.UniformB()]

    def __call__(self):
        return random.choice(self.components)()

    def __str__(self):
        return 'a mixture of branch length distributions'


def get_sorted_eigensystem(M):
    w, vT = np.linalg.eigh(M)
    w_v_pairs = zip(w, vT.T.tolist())
    return list(sorted(w_v_pairs))

def get_stability(D):
    """
    The stability is defined as a bound on norms of perturbation matrices.
    If D is perturbed by a matrix whose Frobenius norm
    is less than the stability, then the spectral split remains unchanged.
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

def process(ntaxa, nseconds, branch_length_sampler):
    """
    The sampling functor returns a branch length and has a string cast.
    @param ntaxa: the number of taxa in the sampled trees
    @param nseconds: allow this many seconds to run or None to run forever
    @param branch_length_sampler: a sampling functor
    @return: a multi-line string that summarizes the results
    """
    data_rows = []
    start_time = time.time()
    try:
        while True:
            elapsed_time = time.time() - start_time
            if nseconds and elapsed_time > nseconds:
                break
            # sample the tree
            tree = TreeSampler.sample_agglomerated_tree(ntaxa)
            # get the atteson bound
            for branch in tree.get_branches():
                branch.length = branch_length_sampler()
            atteson_bound = 0.5 * min(b.length for b in tree.get_branches())
            # get the spectral bound
            D = np.array(tree.get_distance_matrix())
            k = len(D)
            spectral_bound = get_stability(D) / k
            # store the row
            row = [atteson_bound, spectral_bound, tree.get_newick_string()]
            data_rows.append(row)
    except KeyboardInterrupt, e:
        pass
    # make the response
    out = StringIO()
    print >> out, '#', elapsed_time, 'seconds of run time'
    print >> out, '#', ntaxa, 'taxa per tree'
    print >> out, '#', branch_length_sampler
    print >> out, '\t'.join(['atteson.bound', 'spectral.bound', 'newick'])
    for i, row in enumerate(data_rows):
        atteson, spectral, newick = row
        print >> out, '\t'.join([
            str(i+1), str(atteson), str(spectral), '"' + newick + '"'])
    return out.getvalue().strip()

def get_response_content(fs):
    nseconds = 2
    ntaxa = fs.ntaxa
    branch_length_sampler = BranchLengthComboPack()
    return process(ntaxa, nseconds, branch_length_sampler) + '\n'

def main(options):
    assert 0 <= options.nseconds
    assert 4 <= options.ntaxa <= 20
    branch_length_sampler = BranchLengthComboPack()
    print process(options.ntaxa, options.nseconds, branch_length_sampler)

#FIXME use argparse with types instead of using assert

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('--nseconds', dest='nseconds', type='int', default=0,
            help='seconds to run or 0 to run until ctrl-c')
    parser.add_option('--ntaxa', dest='ntaxa', type='int', default=20,
            help='number of taxa in each sampled tree topology')
    options, args = parser.parse_args()
    main(options)
