"""Characterize the efficiency of a distance matrix sampler.

Each distance matrix is sampled in three steps.
First, a nucleotide alignment is sampled from the distribution implied by the tree using a Jukes-Cantor model.
Second, the maximum likelihood distance between each sequence pair is calculated.
Third, the sampled matrix may be rejected if it has elements that are zero or infinity.
"""

import StringIO
import time

import profile

from SnippetUtil import HandlingError
import NewickIO
import FelTree
import DistanceMatrixSampler
import Form

def get_form():
    """
    @return: the body of a form
    """
    # define the default tree string
    tree_string = '(((a:0.05, b:0.05):0.15, c:0.2):0.8, x:1.0, (((m:0.05, n:0.05):0.15, p:0.2):0.8, y:1.0):1.0);'
    tree = NewickIO.parse(tree_string, FelTree.NewickTree)
    formatted_tree_string = NewickIO.get_narrow_newick_string(tree, 60)
    # define the form objects
    form_objects = [
            Form.MultiLine('tree', 'newick tree', formatted_tree_string),
            Form.Integer('length', 'use sequences that are this long', 100, low=1),
            Form.RadioGroup('assumption', 'distance matrix sampling assumption', [
                Form.RadioItem('infinite_alleles', 'infinite alleles', True),
                Form.RadioItem('jukes_cantor', 'Jukes-Cantor model (4 alleles)')]),
            Form.RadioGroup('infinity', 'options for handling distance estimates of infinity', [
                Form.RadioItem('reject_infinity', 'reject matrices with infinite distances'),
                Form.RadioItem('replace_infinity', 'use the value 20 instead of an infinite distance', True)]),
            Form.RadioGroup('zero', 'options for handling distance estimates of zero', [
                Form.RadioItem('reject_zero', 'reject matrices with distances of zero'),
                Form.RadioItem('replace_zero', 'use the value .00001 instead of a distance of zero'),
                Form.RadioItem('remain_zero', 'estimates of zero are unmodified', True)])]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # read the tree
    tree = NewickIO.parse(fs.tree, FelTree.NewickTree)
    # get arbitrarily ordered leaf names
    ordered_names = list(node.name for node in tree.gen_tips())
    # define the distance matrix sampler
    if fs.infinite_alleles:
        sampler = DistanceMatrixSampler.InfiniteAllelesSampler(tree, ordered_names, fs.length)
    elif fs.jukes_cantor:
        sampler = DistanceMatrixSampler.DistanceMatrixSampler(tree, ordered_names, fs.length)
    if fs.reject_infinity:
        sampler.set_inf_replacement(None)
    elif fs.replace_infinity:
        sampler.set_inf_replacement(20)
    if fs.reject_zero:
        sampler.set_zero_replacement(None)
    elif fs.replace_zero:
        sampler.set_zero_replacement(0.00001)
    elif fs.remain_zero:
        sampler.set_zero_replacement(0.0)
    # define the amount of time allotted to the sampler
    allocated_seconds = 2
    # do some sampling, saving a summary but discarding the samples
    start_time = time.clock()
    run_seconds = 0
    for result in sampler.gen_samples_or_none():
        run_seconds = time.clock() - start_time
        if run_seconds > allocated_seconds:
            break
    # define the response
    out = StringIO.StringIO()
    print >> out, 'these are the results for a', run_seconds, 'second run:'
    print >> out, sampler.proposed, 'samples were proposed'
    print >> out, sampler.accepted, 'samples were accepted'
    print >> out, sampler.proposals_with_zero, 'proposals had a distance estimate of zero'
    print >> out, sampler.proposals_with_inf, 'proposals had a distance estimate of infinity'
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()

def main():
    # use the default sequence length
    sequence_length = 100
    # use the default tree
    tree_string = '(((a:0.05, b:0.05):0.15, c:0.2):0.8, x:1.0, (((m:0.05, n:0.05):0.15, p:0.2):0.8, y:1.0):1.0);'
    tree = NewickIO.parse(tree_string, FelTree.NewickTree)
    # get arbitrarily ordered leaf names
    ordered_names = list(node.name for node in tree.gen_tips())
    # create the sampler
    sampler = DistanceMatrixSampler.InfiniteAllelesSampler(tree, ordered_names, sequence_length)
    sampler.set_inf_replacement(20)
    sampler.set_zero_replacement(0.0)
    # do some sampling, saving a summary but discarding the samples
    allocated_seconds = 2
    start_time = time.clock()
    run_seconds = 0
    for result in sampler.gen_samples_or_none():
        run_seconds = time.clock() - start_time
        if run_seconds > allocated_seconds:
            break
    # define the response
    out = StringIO.StringIO()
    print 'these are the results for a', run_seconds, 'second run:'
    print sampler.proposed, 'samples were proposed'
    print sampler.accepted, 'samples were accepted'
    print sampler.proposals_with_zero, 'proposals had a distance estimate of zero'
    print sampler.proposals_with_inf, 'proposals had a distance estimate of infinity'

if __name__ == '__main__':
    profile.run('main()')

