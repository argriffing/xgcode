"""Characterize the efficiency of a tree reconstruction sampler.

Each distance matrix is sampled in three steps.
First, a nucleotide alignment is sampled from the distribution
implied by the tree using a Jukes-Cantor model.
Second, the maximum likelihood distance
between each sequence pair is calculated.
Third, the sampled matrix may be rejected
if it has elements that are zero or infinity.
"""

from StringIO import StringIO
import time

from SnippetUtil import HandlingError
import SnippetUtil
import NewickIO
import FelTree
import DistanceMatrixSampler
import NeighborhoodJoining
import Clustering
from Form import RadioItem
import Form
import FormOut

def get_form():
    """
    @return: the body of a form
    """
    # define the tree string
    tree_string = '(((a:0.05, b:0.05):0.15, c:0.2):0.8, x:1.0, (((m:0.05, n:0.05):0.15, p:0.2):0.8, y:1.0):1.0);'
    tree = NewickIO.parse(tree_string, FelTree.NewickTree)
    formatted_tree_string = NewickIO.get_narrow_newick_string(tree, 60)
    # define the object list
    form_objects = [
            Form.MultiLine('tree', 'tree',
                formatted_tree_string),
            Form.Integer('sequence_length', 'use sequences that are this long',
                100, low=1),
            Form.RadioGroup('assumption', 'distance matrix sampling model', [
                RadioItem('infinite_alleles', 'infinite alleles', True),
                RadioItem('jukes_cantor', 'Jukes-Cantor model (4 alleles)')]),
            Form.RadioGroup('infinity', 'matrices with infinite distances', [
                RadioItem('reject_infinity', 'reject these matrices', True),
                RadioItem('replace_infinity', 'use 20 instead')]),
            Form.RadioGroup('zero', 'matrices with zero distances', [
                RadioItem('reject_zero', 'reject these matrices'),
                RadioItem('replace_zero', 'use .00001 instead'),
                RadioItem('remain_zero', 'use 0 unmodified', True)]),
            Form.RadioGroup('criterion', 'tree reconstruction criterion', [
                RadioItem('sign', 'spectral sign approximation', True),
                RadioItem('nj', 'neighbor joining'),
                RadioItem('random', 'random bipartition')])]
    # return the object list
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # read the tree
    tree = NewickIO.parse(fs.tree, FelTree.NewickTree)
    # read the sequence length
    sequence_length = fs.sequence_length
    # get arbitrarily ordered leaf names
    ordered_names = list(node.name for node in tree.gen_tips())
    # read the criterion string, creating the splitter object
    if fs.sign:
        splitter = Clustering.StoneSpectralSignDMS()
    elif fs.nj:
        splitter = Clustering.NeighborJoiningDMS()
    elif fs.random:
        splitter = Clustering.RandomDMS()
    # define the distance matrix sampler
    if fs.infinite_alleles:
        sampler = DistanceMatrixSampler.InfiniteAllelesSampler(
                tree, ordered_names, sequence_length)
    elif fs.jukes_cantor:
        sampler = DistanceMatrixSampler.DistanceMatrixSampler(
                tree, ordered_names, sequence_length)
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
    allocated_seconds = 1
    # get distance matrices until we run out of time
    distance_matrices = []
    start_time = time.clock()
    sampling_seconds = 0
    for result in sampler.gen_samples_or_none():
        # if the result was accepted then add the distance matrix
        if result is not None:
            sequence_list, D = result
            distance_matrices.append(D)
        # see if we need to stop sampling
        sampling_seconds = time.clock() - start_time
        if sampling_seconds >= allocated_seconds:
            break
    # reconstruct trees until we run out of time
    start_time = time.clock()
    reconstructing_seconds = 0
    reconstructed_tree_count = 0
    for D in distance_matrices:
        # reconstruct a tree using the method of choice
        tree_builder = NeighborhoodJoining.TreeBuilder(
                D, ordered_names, splitter)
        tree_builder.set_fallback_name('nj')
        try:
            query_tree = tree_builder.build()
        except NeighborhoodJoining.NeighborhoodJoiningError, e:
            raise HandlingError(e)
        reconstructed_tree_count += 1
        # see if we need to stop reconstructing the trees
        reconstructing_seconds = time.clock() - start_time
        if reconstructing_seconds >= allocated_seconds:
            break
    # define the response
    out = StringIO()
    if distance_matrices:
        print >> out, 'seconds to sample', len(distance_matrices),
        print >> out, 'distance matrices:', sampling_seconds
        if reconstructed_tree_count:
            print >> out, 'seconds to reconstruct', reconstructed_tree_count,
            print >> out, 'trees:', reconstructing_seconds
        else:
            print >> out, 'no trees could be reconstructed',
            print >> out, 'in a reasonable amount of time'
    else:
        print >> out, 'no distance matrices could be sampled'
        print >> out, 'in a reasonable amount of time'
        print >> out, sampler.proposed,
        print >> out, 'distance matrices were proposed but were rejected'
        print >> out, sampler.proposals_with_zero,
        print >> out, 'proposed distance matrices had estimates of zero'
        print >> out, sampler.proposals_with_inf,
        print >> out, 'proposed distance matrices had estimates of infinity'
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()

