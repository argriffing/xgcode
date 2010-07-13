"""Evaluate the quality of various distance matrix tree reconstruction methods.

This evaluation procedure is somewhat complicated.
Several nucleotide alignments are simulated
using the JC69 model and the original tree.
Distance matrices are estimated from these alignments using maximum likelihood,
and distance matrices containing an element that is infinite are rejected.
The tree reconstruction method then builds a tree
from each estimated distance matrix.
The reconstructed trees are then compared to the original tree,
and the distribution of the number of implied partition errors is reported.
The exact bipartition criterion is a matrix function by Eric Stone.
"""

from StringIO import StringIO

from SnippetUtil import HandlingError
import MatrixUtil
import NewickIO
import FelTree
import Clustering
import NeighborhoodJoining
import JC69
import PhyLikelihood
import RateMatrix
from Form import RadioItem
from Form import CheckItem
import Form
import FormOut

#FIXME use distance matrix sampler module
#FIXME use const data

class DistanceMatrixSampler:
    """
    Sample estimated distance matrices,
    rejecting those with infinite branch lengths.
    The complexity of generating the samples is also estimated,
    which is important because otherwise the rejection sampling
    could get stuck in a loop if every sample is rejected.
    """

    def __init__(self, tree, ordered_names, sequence_length):
        assert len(list(tree.gen_tips())) == len(ordered_names)
        self.tree = tree
        self.ordered_names = ordered_names
        self.sequence_length = sequence_length
        self.requested_matrix_count = 0
        self.accepted_sample_count = 0
        # Initialize the number of samples rejected
        # because of an infinitely long branch length estimate.
        self.rejected_inf_sample_count = 0
        # Initialize the number of samples rejected
        # because of a branch length estimate of zero.
        self.rejected_zero_sample_count = 0

    def get_sampling_error_message(self):
        accepted = self.accepted_sample_count
        rejected_inf = self.rejected_inf_sample_count
        rejected_zero = self.rejected_zero_sample_count
        rejected = rejected_inf + rejected_zero
        total = accepted + rejected
        msg_c = None
        if total:
            msg_c = '%d of %d samples accepted' % (accepted, total)
            if accepted < rejected_inf:
                msg_c += '; use shorter branch lengths or longer sequences'
            elif accepted < rejected_zero:
                msg_c += '; use longer branch lengths or longer sequences'
        msg_a = 'the distance matrix sampling procedure takes too long '
        msg_b = 'for these settings'
        error_message = msg_a + msg_b
        if msg_c:
            error_message += ' (%s)' % msg_c
        return error_message

    def get_rejected_sample_count(self):
        return self.rejected_inf_sample_count + self.rejected_zero_sample_count

    def get_acceptance_probability(self):
        """
        This is for progress bar stuff.
        @return: an optimistic acceptance probability for rejection sampling
        """
        total_samples = (
                self.accepted_sample_count + self.get_rejected_sample_count())
        if total_samples < 100:
            # if not enough samples have been taken then be optimistic
            return 1.0
        else:
            # If a reasonable number of samples have been taken
            # then be realistic.
            return self.accepted_sample_count / float(total_samples)

    def get_complexity(self):
        """
        This is for progress bar stuff.
        Return the predicted total number of steps required,
        for some step granularity.
        @return: the predicted total number of steps
        """
        # if all of the samples are rejected then the complexity is infinite
        acceptance_probability = self.get_acceptance_probability()
        if not acceptance_probability:
            return float('inf')
        # If there is some predicted probability of accepting a sample
        # then make a guess.
        n = len(self.ordered_names)
        steps_per_sample = n * n * self.sequence_length
        required_accepted_samples = self.requested_matrix_count
        samples_per_accepted_sample = 1.0 / acceptance_probability
        # compute the number of steps
        steps = steps_per_sample
        steps *= samples_per_accepted_sample * required_accepted_samples
        return steps

    def gen_distance_matrices(self, count, max_steps):
        """
        Yield (ordered sequence list, distance matrix) pairs .
        The generator will stop if it sees that it cannot meet its goal
        in the allotted number of steps.
        @param count: the requested number of distance matrices
        @param max_steps: an upper bound on the allowed number of steps
        """
        # define the jukes cantor rate matrix
        dictionary_rate_matrix = RateMatrix.get_jukes_cantor_rate_matrix()
        ordered_states = list('ACGT')
        row_major_rate_matrix = MatrixUtil.dict_to_row_major(
                dictionary_rate_matrix, ordered_states, ordered_states)
        model = RateMatrix.RateMatrix(row_major_rate_matrix, ordered_states)
        # record the requested number of samples
        self.requested_matrix_count = count
        # do some rejection sampling
        while True:
            if self.get_complexity() >= max_steps:
                break
            if self.accepted_sample_count >= count:
                break
            # simulate an alignment from the tree
            alignment = PhyLikelihood.simulate_alignment(
                    self.tree, model, self.sequence_length)
            # extract the ordered list of sequences from the alignment object
            name_to_sequence = dict(zip(alignment.headers, alignment.sequences))
            sequence_list = [name_to_sequence[name]
                    for name in self.ordered_names]
            # get the estimated distance matrix
            distance_matrix = JC69.get_ML_distance_matrix(sequence_list)
            # look for degeneracies
            has_zero_off_diagonal = False
            has_inf_off_diagonal = False
            for i, row in enumerate(distance_matrix):
                for j, value in enumerate(row):
                    if i != j:
                        if value == 0.0:
                            has_zero_off_diagonal = True
                        if value == float('inf'):
                            has_inf_off_diagonal = True
            if has_zero_off_diagonal:
                self.rejected_zero_sample_count += 1
            elif has_inf_off_diagonal:
                self.rejected_inf_sample_count += 1
            else:
                self.accepted_sample_count += 1
                yield sequence_list, distance_matrix


def get_default_original_tree():
    tree_string = '(((a:0.05, b:0.05):0.15, c:0.2):0.8, x:1.0, (((m:0.05, n:0.05):0.15, p:0.2):0.8, y:1.0):1.0);'
    tree = NewickIO.parse(tree_string, FelTree.NewickTree)
    for node in tree.preorder():
        blen = node.get_branch_length()
        if blen is not None:
            node.set_branch_length(blen * 0.5)
    return tree

def get_form():
    """
    @return: the body of a form
    """
    # define the default tree string
    tree = get_default_original_tree()
    formatted_tree_string = NewickIO.get_narrow_newick_string(tree, 60)
    # define the form objects
    form_objects = [
            Form.MultiLine('tree', 'original tree with branch lengths',
                formatted_tree_string),
            Form.Integer('iterations', 'reconstruct this many trees',
                10, low=1),
            Form.Integer('length', 'use sequences that are this long',
                100, low=2),
            Form.RadioGroup('criterion', 'bipartition function', [
                RadioItem('exact', 'exact criterion'),
                RadioItem('sign', 'spectral sign approximation', True),
                RadioItem('threshold', 'spectral threshold approximation'),
                RadioItem('nj', 'neighbor joining criterion'),
                RadioItem('random', 'random bipartition')]),
            Form.RadioGroup('recourse', 'recourse for degenerate partitions', [
                RadioItem('njrecourse', 'neighbor joining', True),
                RadioItem('halvingrecourse', 'leaf stem length halving')]),
            Form.CheckGroup('output_options', 'extra output option', [
                CheckItem('showtrees', 'show reconstructed tree topologies')])]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # read the criterion string, creating the splitter object
    if fs.exact:
        splitter = Clustering.StoneExactDMS()
    elif fs.sign:
        splitter = Clustering.StoneSpectralSignDMS()
    elif fs.threshold:
        splitter = Clustering.StoneSpectralThresholdDMS()
    elif fs.nj:
        splitter = Clustering.NeighborJoiningDMS()
    elif fs.random:
        splitter = Clustering.RandomDMS()
    # read the original tree
    tree = NewickIO.parse(fs.tree, FelTree.NewickTree)
    # Make sure that the splitter object is appropriate for the number
    # of taxa and the number of tree reconstructions.
    ntaxa = len(list(tree.gen_tips()))
    if splitter.get_complexity(ntaxa) * fs.iterations > 1000000:
        msg_a = 'use a faster bipartition function, fewer taxa, '
        msg_b = 'or fewer tree reconstructions'
        raise HandlingError(msg_a + msg_b)
    # sample a bunch of sequences
    ordered_names = [node.name for node in tree.gen_tips()]
    sampler = DistanceMatrixSampler(tree, ordered_names, fs.length)
    # simulate a bunch of distance matrices and reconstruct the trees
    mismatch_count_tree_pairs = []
    error_count_histogram = {}
    max_steps = 1000000
    for sequence_list, distance_matrix in sampler.gen_distance_matrices(
            fs.iterations, max_steps):
        # create the tree builder
        tree_builder = NeighborhoodJoining.ValidatingTreeBuilder(
                distance_matrix, ordered_names, splitter)
        # Read the recourse string and set the corresponding method
        # in the tree builder.
        if fs.njrecourse:
            tree_builder.set_fallback_name('nj')
        elif fs.halvingrecourse:
            tree_builder.set_fallback_name('halving')
        # set parameters of the tree validating tree builder
        tree_builder.set_original_tree(tree)
        # build the tree
        reconstructed_tree = tree_builder.build()
        # note the number of partition errors during the reconstruction
        mismatch_count = tree_builder.get_mismatch_count()
        if mismatch_count not in error_count_histogram:
            error_count_histogram[mismatch_count] = 0
        error_count_histogram[mismatch_count] += 1
        # If we are saving the reconstructed trees
        # then remove branch lengths and add to the tree list.
        if fs.showtrees:
            for node in reconstructed_tree.preorder():
                node.set_branch_length(None)
            mismatch_count_tree_pair = (mismatch_count, reconstructed_tree)
            mismatch_count_tree_pairs.append(mismatch_count_tree_pair)
    # See if we bailed early because
    # the sampling was predicted to take too long.
    if sampler.accepted_sample_count < fs.iterations:
        raise HandlingError(sampler.get_sampling_error_message())
    # define the response
    out = StringIO()
    print >> out, 'partition error count frequencies:'
    max_mismatch_count = max(error_count_histogram)
    for i in range(max_mismatch_count + 1):
        frequency = error_count_histogram.get(i, 0)
        print >> out, i, ':', frequency
    if fs.showtrees:
        print >> out, ''
        print >> out, 'reconstructed tree topologies with mismatch counts:'
        for mismatch_count, tree in sorted(mismatch_count_tree_pairs):
            print >> out, NewickIO.get_newick_string(tree), mismatch_count
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()


class Simulation:
    """
    This class represents a simulation run of a reconstruction method.
    It is meant to be used when run from the command line.
    """

    def __init__(self, splitter, fallback_name, description):
        """
        The description is of the method used to split the distance matrix.
        @param splitter: a distance matrix splitter
        @param fallback_name: the name
        @param description: a description of the method 
        """
        # These simulation parameters are set
        # at initialization time.
        self.splitter = splitter
        self.fallback_name = fallback_name
        self.description = description
        # These simulation parameters are set
        # after the object has been initialized.
        self.sequence_length = None
        self.step_limit = None
        self.original_tree = None
        self.reconstruction_count = None
        # this is internal state data
        self.histogram = {}

    def add_error_count(self, error_count):
        """
        Add an error count.
        The error count is the number of partition errors
        in a reconstructed tree relative to the original tree.
        @param error_count: the number of partition errors
        """
        if error_count not in self.histogram:
            self.histogram[error_count] = 0
        self.histogram[error_count] += 1

    def get_count_list(self):
        """
        Get a list of error counts.
        The first element of the returned list
        is the number of times that no errors occurred.
        The second element
        is the number of times that one error occurred.
        The length of the list is equal to
        the number of errors in the reconstruction with the most errors.
        @return: a list of error counts
        """
        max_error_count = max(self.histogram)
        return [self.histogram.get(i, 0) for i in range(max_error_count + 1)]

    def get_histogram_string(self):
        """
        Return a multi-line string.
        It summarizes the quality of the trees
        reconstructed during the simulation
        @return: a multi-line string
        """
        out = StringIO()
        for i, count in enumerate(self.get_count_list()):
            print >> out, i, ':', count
        return out.getvalue().strip()

    def set_original_tree(self, original_tree):
        """
        @param original_tree: the true tree with branch lengths
        """
        self.original_tree = original_tree

    def set_reconstruction_count(self, reconstruction_count):
        """
        @param reconstruction_count: the number of simulations to do
        """
        self.reconstruction_count = reconstruction_count

    def set_step_limit(self, step_limit):
        """
        @param step_limit: a cap on the number of steps allowed
        """
        self.step_limit = step_limit

    def set_sequence_length(self, sequence_length):
        """
        @param sequence_length: this is the length of the simulated sequences
        """
        self.sequence_length = sequence_length

    def run(self):
        # simulate a bunch of distance matrices
        ordered_names = [node.name for node in self.original_tree.gen_tips()]
        sampler = DistanceMatrixSampler(
                self.original_tree, ordered_names, self.sequence_length)
        for sequence_list, distance_matrix in sampler.gen_distance_matrices(
                self.reconstruction_count, self.step_limit):
            # create the tree builder
            tree_builder = NeighborhoodJoining.ValidatingTreeBuilder(
                    distance_matrix, ordered_names, self.splitter)
            # set parameters of the tree validating tree builder
            tree_builder.set_fallback_name(self.fallback_name)
            tree_builder.set_original_tree(self.original_tree)
            # build the tree
            try:
                reconstructed_tree = tree_builder.build()
            except NeighborhoodJoining.NeighborhoodJoiningError, e:
                print 'neighborhood joining error:', e
                print 'simulated sequence list:'
                for sequence in sequence_list:
                    print sequence
            # note the number of partition errors during the reconstruction
            self.add_error_count(tree_builder.get_mismatch_count())
        # See if we bailed early
        # because the sampling was predicted to take too long.
        if sampler.accepted_sample_count < self.reconstruction_count:
            raise HandlingError(sampler.get_sampling_error_message())


def main():
    """
    Run some tree reconstructions from the command line.
    """
    # initialize the simulation objects
    sims = [
        Simulation(Clustering.NeighborJoiningDMS(),
            'nj', 'neighbor joining'),
        Simulation(Clustering.RandomDMS(),
            'nj', 'random partitioning'),
        Simulation(Clustering.StoneExactDMS(),
            'nj', 'exact criterion with neighbor joining fallback'),
        #Simulation(Clustering.StoneExactDMS(),
        #'halving', 'exact criterion with stem halving fallback'),
        Simulation(Clustering.StoneSpectralSignDMS(),
            'nj', 'spectral sign cut with neighbor joining fallback')
        #Simulation(Clustering.StoneSpectralSignDMS(),
        #'halving', 'spectral sign cut with stem halving fallback')
        ]
    # define the simulation parameters
    tree = get_default_original_tree()
    reconstruction_count = 1000
    sequence_length = 100
    step_limit_per_method = 10000000
    # set the simulation parameters
    for sim in sims:
        sim.set_original_tree(get_default_original_tree())
        sim.set_reconstruction_count(reconstruction_count)
        sim.set_step_limit(step_limit_per_method)
        sim.set_sequence_length(sequence_length)
    # show the simulation parameters
    print 'simulation parameters:'
    print 'original tree:', NewickIO.get_newick_string(tree)
    print 'reconstruction count:', reconstruction_count
    print 'sequence length:', sequence_length
    # run the simulations
    print 'running the simulations...'
    for sim in sims:
        print 'running "%s"...' % sim.description
        try:
            sim.run()
        except HandlingError, e:
            print 'Error:', e
    # print the simulation data
    print 'simulation results:'
    for sim in sims:
        print sim.description + ':'
        print sim.get_histogram_string()

if __name__ == '__main__':
    main()
