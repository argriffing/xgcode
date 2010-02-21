"""Evaluate distance matrix tree inference, giving more weight to deep splits.

Evaluate distance matrix tree building methods,
giving more weight to deep splits.
This evaluation procedure is somewhat complicated.
Several nucleotide alignments are simulated using the JC69 model
and the original tree.
Distance matrices are estimated from these alignments using maximum likelihood,
and distance matrices containing an element that is infinite are rejected.
The tree reconstruction method then builds a tree
from each estimated distance matrix.
The reconstructed trees are then compared to the original tree,
and the weight of the partition errors is reported.
"""

from StringIO import StringIO
import optparse
import time
import math

from SnippetUtil import HandlingError
import SnippetUtil
import Util
import NewickIO
import FelTree
import Clustering
import NeighborhoodJoining
import DistanceMatrixSampler
import JC69
import HtmlTable
import TreeComparison
import RUtil
import Form
import Progress


class OptionError(Exception):
    pass


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
            Form.Integer('iterations', 'reconstruct this many trees', 10),
            Form.Integer('length', 'use sequences that are this long', 100),
            Form.RadioGroup('criterion', 'bipartition function', [
                Form.RadioItem('exact', 'exact criterion'),
                Form.RadioItem('sign', 'spectral sign approximation', True),
                Form.RadioItem('nj', 'neighbor joining criterion'),
                Form.RadioItem('random', 'random bipartition')])]
    return form_objects

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
    elif fs.nj:
        splitter = Clustering.NeighborJoiningDMS()
    elif fs.random:
        splitter = Clustering.RandomDMS()
    # read the original tree
    tree = NewickIO.parse(fs.tree, FelTree.NewickTree)
    # define the maximum number of steps we want
    max_steps = 1000000
    # Make sure that the splitter object is appropriate
    # for the number of taxa and the number of tree reconstructions.
    ntaxa = len(list(tree.gen_tips()))
    if splitter.get_complexity(ntaxa) * fs.iterations > max_steps:
        msg_a = 'use a faster bipartition function, '
        msg_b = 'fewer taxa, or fewer tree reconstructions'
        raise HandlingError(msg_a + msg_b)
    # define the simulation parameters
    sim = Simulation(splitter, 'nj', 'cgi tree building simulation')
    sim.set_original_tree(tree)
    sim.set_step_limit(max_steps)
    # define an arbitrary but consistent ordering of the taxa
    ordered_names = [node.name for node in tree.gen_tips()]
    # attempt to simulate a bunch of distance matrices
    sampler = DistanceMatrixSampler.DistanceMatrixSampler(
            tree, ordered_names, fs.length)
    distance_matrices = []
    for result in sampler.gen_samples_or_none():
        # if a proposal was accepted then add it to the list
        if result:
            sequence_list, distance_matrix = result
            distance_matrices.append(distance_matrix)
        # if enough accepted samples have been generated then stop sampling
        remaining_acceptances = fs.iterations - len(distance_matrices)
        if not remaining_acceptances:
            break
        # If the remaining number of computrons is predicted
        # to be too much then stop.
        if sampler.get_remaining_computrons(remaining_acceptances) > max_steps:
            msg_a = 'this combination of parameters '
            msg_b = 'is predicted to take too long'
            raise HandlingError(msg)
    sim.run(distance_matrices, ordered_names)
    # define the response
    out = StringIO()
    print >> out, 'partition error count frequencies:'
    print >> out, sim.get_histogram_string()
    print >> out, ''
    print >> out, 'weighted partition errors:', sim.get_deep_loss()
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
        The description is of the distance matrix splitting method.
        @param splitter: a distance matrix splitter
        @param fallback_name: the name
        @param description: a description of the splitting method
        """
        # these simulation parameters are set at initialization time
        self.splitter = splitter
        self.fallback_name = fallback_name
        self.description = description
        # These simulation parameters are set
        # after the object has been initialized.
        self.verbose = False
        self.step_limit = None
        self.original_tree = None
        # For each reconstructed tree,
        # record two distances to the original tree.
        self.error_counts = []
        self.loss_values = []
        self.max_error_counts = []
        self.max_loss_values = []
        # timing information
        self.start_time = None
        self.stop_time = None

    def get_uniform_loss(self):
        """
        Return the averaged proportion of partition errors.
        It is averaged over all reconstructed trees.
        @return: averaged proportion of partition errors
        """
        numerator = sum(self.error_counts)
        denominator = sum(self.max_error_counts)
        if denominator:
            if numerator > denominator:
                raise HandlingError('uniform loss normalization error')
            return numerator / float(denominator)
        else:
            return float('inf')

    def get_deep_loss(self):
        """
        Return the averaged proportion of partition errors.
        It is averaged over all reconstructed trees.
        @return: averaged proportion of partition errors
        """
        numerator = sum(self.loss_values)
        denominator = sum(self.max_loss_values)
        if denominator:
            if numerator > denominator:
                raise HandlingError('deep loss normalization error')
            return numerator / float(denominator)
        else:
            return float('inf')

    def get_count_list(self):
        """
        The first element of the returned list
        is the number of times that no errors occurred.
        The second element is the number of times that one error occurred.
        The length of the list is equal to the number of errors
        in the reconstruction with the most errors.
        @return: a list of error counts
        """
        max_error_count = max(self.error_counts)
        count_list = [0] * (max_error_count + 1)
        for count in self.error_counts:
            count_list[count] += 1
        return count_list

    def get_histogram_string(self):
        """
        Get a histogram.
        Return a multi-line string summarizing the quality of the trees
        reconstructed during the simulation.
        @return: multi-line histogram string
        """
        out = StringIO()
        for i, count in enumerate(self.get_count_list()):
            print >> out, i, ':', count
        return out.getvalue().strip()

    def set_verbose(self, verbose=True):
        """
        @param verbose: True when we want verbose output
        """
        self.verbose = verbose

    def set_original_tree(self, original_tree):
        """
        @param original_tree: the true tree with branch lengths
        """
        self.original_tree = original_tree

    def set_step_limit(self, step_limit):
        """
        @param step_limit: limit the number of steps allowed in the computation
        """
        self.step_limit = step_limit

    def get_running_time(self):
        """
        @return: the number of seconds it took to run the simulation
        """
        if self.start_time is None:
            raise HandlingError('the simulation has not been started')
        if self.stop_time is None:
            msg = 'the simulation was not successfully completed'
            raise HandlingError(msg)
        return self.stop_time - self.start_time

    def get_normalized_error_counts(self):
        """
        @return: a list of normalized error counts from the completed run
        """
        normalized_error_counts = []
        for error_count, max_error_count in zip(
                self.error_counts, self.max_error_counts):
            numerator = float(error_count)
            denominator = float(max_error_count)
            assert numerator <= denominator
            normalized_error_counts.append(numerator / denominator)
        return normalized_error_counts

    def get_normalized_loss_values(self):
        """
        @return: a list of normalized loss values from the completed run
        """
        normalized_loss_values = []
        for loss_value, max_loss_value in zip(
                self.loss_values, self.max_loss_values):
            numerator = float(loss_value)
            denominator = float(max_loss_value)
            assert numerator <= denominator
            normalized_loss_values.append(numerator / denominator)
        return normalized_loss_values

    def run(self, distance_matrices, ordered_names):
        """
        This function stores the losses for each reconstruction.
        @param distance_matrices: a sequence of distance matrices
        @param ordered_names: order of taxa in the distance matrix
        """
        if self.start_time is not None:
            msg = 'each simulation object should be run only once'
            raise HandlingError(msg)
        if not distance_matrices:
            raise HandlingErrror('no distance matrices were provided')
        tip_name_set = set(node.name for node in self.original_tree.gen_tips()) 
        if tip_name_set != set(ordered_names):
            raise HandlingError('leaf name mismatch')
        self.start_time = time.time()
        # Define the reference tree and its maximum cost
        # under different loss functions.
        reference_tree = self.original_tree
        max_error_count = TreeComparison.get_nontrivial_split_count(
                reference_tree)
        max_loss_value = TreeComparison.get_weighted_split_count(
                reference_tree)
        for distance_matrix in distance_matrices:
            # create the tree builder
            tree_builder = NeighborhoodJoining.TreeBuilder(
                    distance_matrix, ordered_names, self.splitter)
            # set parameters of the validating tree builder
            tree_builder.set_fallback_name(self.fallback_name)
            # build the tree
            try:
                query_tree = tree_builder.build()
            except NeighborhoodJoining.NeighborhoodJoiningError, e:
                raise HandlingError(e)
            # Note the number and weight of partition errors
            # during the reconstruction.
            error_count = TreeComparison.get_split_distance(
                    query_tree, reference_tree)
            loss_value = TreeComparison.get_weighted_split_distance(
                    query_tree, reference_tree)
            # make sure that the summary is internally consistent
            assert error_count <= max_error_count, (
                    error_count, max_error_count)
            assert loss_value <= max_loss_value, (
                    loss_value, max_loss_value)
            # save the reconstruction characteristics to use later
            self.error_counts.append(error_count)
            self.loss_values.append(loss_value)
            self.max_error_counts.append(max_error_count)
            self.max_loss_values.append(max_loss_value)
        self.stop_time = time.time()

def get_tree_and_remark(options):
    """
    @param options: an object from optparse
    @return: a tree and a string that is a comment about the tree
    """
    if options.inline_tree:
        tree_string = options.inline_tree
        tree = NewickIO.parse(tree_string, FelTree.NewickTree)
        tree_remark = options.inline_tree
    elif options.tree_filename:
        fin = open(options.tree_filename)
        tree_string = fin.read()
        fin.close()
        tree = NewickIO.parse(tree_string, FelTree.NewickTree)
        tree_remark = options.tree_filename
    else:
        tree = get_default_original_tree()
        tree_remark = 'a default tree'
    return (tree, tree_remark)

def R_helper(nj_losses, ss_losses):
    """
    @param nj_losses: a list of neighbor joining loss values
    @param ss_losses: a list of spectral sign loss values
    @return: the contents of an R file as a string
    """
    #nj_log_losses = [math.log(x) for x in nj_losses]
    #ss_log_losses = [math.log(x) for x in ss_losses]
    nj_string = ', '.join(str(x) for x in nj_losses)
    ss_string = ', '.join(str(x) for x in ss_losses)
    arr = [
            'nj <- c(%s)' % nj_string,
            'ss <- c(%s)' % ss_string
            ]
    return '\n'.join(arr)

def do_hard_coded_analysis_a(tree, tree_remark):
    """
    Do a hardcoded analysis of tree reconstruction methods.
    Make a bunch of R files.
    @param tree: a tree object
    @param tree_remark: a string that is a comment about the tree
    """
    # define an arbitrary order for the names of the leaves of the tree
    ordered_names = list(node.name for node in tree.gen_tips())
    # use 1000 replicates
    reconstruction_count = 1000
    # Make R files for reconstruction results
    # from sequences 100 and 500 nucleotides long.
    for sequence_length in (100, 500):
        # sample distance matrices
        print 'sampling', reconstruction_count, 'distance matrices'
        print 'from alignments of length', sequence_length
        sampler = DistanceMatrixSampler.DistanceMatrixSampler(
                tree, ordered_names, sequence_length)
        distance_matrices = []
        for result in sampler.gen_samples_or_none():
            # if the proposal was rejected then try again
            if not result:
                continue
            # add the accepted distance matrix sample to the list
            sequence_list, distance_matrix = result
            distance_matrices.append(distance_matrix)
            # stop when we have generated enough distance matrices
            if len(distance_matrices) == reconstruction_count:
                break
        # run both neighbor joining and spectral sign clustering
        sims = [
                Simulation(Clustering.NeighborJoiningDMS(),
                    'nj', 'neighbor joining'),
                Simulation(Clustering.StoneSpectralSignDMS(),
                    'nj', 'spectral sign')]
        for sim in sims:
            print 'reconstructing', len(distance_matrices), 'trees'
            print 'using', sim.description
            sim.set_original_tree(tree)
            sim.run(distance_matrices, ordered_names)
        # consider the neighbor joining and the spectral sign results
        nj_sim, ss_sim = sims
        # write the uniform loss function comparison R script
        script_contents = R_helper(nj_sim.get_normalized_error_counts(),
                ss_sim.get_normalized_error_counts())
        filename = 'uniform_%d.R' % sequence_length
        fout = open(filename, 'w')
        print >> fout, script_contents
        fout.close()
        # write the weighted loss function comparison R script
        script_contents = R_helper(nj_sim.get_normalized_loss_values(),
                ss_sim.get_normalized_loss_values())
        filename = 'weighted_%d.R' % sequence_length
        fout = open(filename, 'w')
        print >> fout, script_contents
        fout.close()

def do_hard_coded_analysis_b(tree, tree_remark):
    """
    Do a hardcoded analysis of tree reconstruction methods.
    Make R files of ordered reconstruction losses.
    @param tree: a tree object
    @param tree_remark: a string that is a comment about the tree
    """
    # define an arbitrary order for the names of the leaves of the tree
    ordered_names = list(node.name for node in tree.gen_tips())
    # use some replicates
    reconstruction_count = 100
    # Make R files for reconstruction results from sequences
    # of some number of nucleotides in length.
    sequence_length = 2000
    # define the tree reconstruction methods to be used
    sims = [
            Simulation(Clustering.NeighborJoiningDMS(),
                'nj', 'neighbor joining'),
            Simulation(Clustering.StoneSpectralSignDMS(),
                'nj', 'spectral sign')]
    # set tree reconstruction parameters
    for sim in sims:
        sim.set_original_tree(tree)
    # initialize the distance matrix sampler
    sampler = DistanceMatrixSampler.InfiniteAllelesSampler(
            tree, ordered_names, sequence_length)
    sampler.set_inf_replacement(20.0)
    sampler.set_zero_replacement(0.0)
    # start the progress bar
    pbar = Progress.Bar(1.0)
    # sample some distance matrices
    distance_matrix_start_time = time.time()
    distance_matrices = []
    for result in sampler.gen_samples_or_none():
        # if we got a result then update the distance matrix list
        if result:
            sequence_list, D = result
            distance_matrices.append(D)
        # Update the progressbar regardless of whether or not
        # the proposal was accepted.
        remaining_acceptances = reconstruction_count - len(distance_matrices)
        numerator = sampler.get_completed_proposals()
        denominator = numerator + sampler.get_remaining_proposals(
                remaining_acceptances)
        dms_fraction = float(numerator) / float(denominator)
        dms_total = 1.0 / (1 + len(sims))
        pbar.update(dms_fraction * dms_total)
        # if we have enough samples then break the loop
        if not remaining_acceptances:
            break
    distance_matrix_seconds = time.time() - distance_matrix_start_time
    # reconstruct trees using various methods
    reconstruction_seconds = []
    for i, sim in enumerate(sims):
        reconstruction_start_time = time.time()
        print 'reconstructing', len(distance_matrices), 'trees'
        print 'using', sim.description
        sim.run(distance_matrices, ordered_names)
        pbar.update(float(i+2) / float(1 + len(sims)))
        reconstruction_seconds.append(time.time() - reconstruction_start_time)
    # stop the progress bar
    pbar.finish()
    # consider the neighbor joining and the spectral sign results
    nj_sim, ss_sim = sims
    # extract the simulation data
    label_list_pairs = [
            ('nj.unweighted', nj_sim.get_normalized_error_counts()),
            ('ss.unweighted', ss_sim.get_normalized_error_counts()),
            ('nj.weighted', nj_sim.get_normalized_loss_values()),
            ('ss.weighted', ss_sim.get_normalized_loss_values())]
    labels, transposed_table = zip(*label_list_pairs)
    table = zip(*transposed_table)
    table_string = RUtil.get_table_string(table, labels)
    # write the table
    filename = 'out3.table'
    fout = open(filename, 'w')
    print >> fout, '# tree source:', tree_remark
    print >> fout, '# number of taxa:', len(ordered_names)
    print >> fout, '# sampled distance matrices:', len(distance_matrices)
    print >> fout, '# seconds elapsed for sampling:', distance_matrix_seconds
    print >> fout, '# sites per sequence:', sequence_length
    for sim, seconds in zip(sims, reconstruction_seconds):
        msg_a = '# seconds elapsed for tree reconstruction using '
        msg_b = sim.description + ': ' + str(seconds)
        print >> fout, msg_a  + msg_b
    print >> fout, table_string
    fout.close()
    print 'wrote', filename

def do_command_line_analysis(options):
    """
    Print some stuff to stdout, and show a progress bar on stderr.
    @param options: an object from optparse
    """
    # load the tree, using the default tree if no filename was provided
    tree, tree_remark = get_tree_and_remark(options)
    # initialize the simulation objects
    sims = [
        Simulation(Clustering.NeighborJoiningDMS(),
            'nj', 'neighbor joining'),
        Simulation(Clustering.StoneSpectralSignDMS(),
            'nj', 'spectral sign cut with neighbor joining fallback'),
        Simulation(Clustering.RandomDMS(),
            'nj', 'random partitioning')]
    # possibly add the slow simulation
    if options.use_exact:
        sims.append(Simulation(Clustering.StoneExactDMS(),
            'nj', 'exact criterion with neighbor joining fallback'))
    # define the simulation parameters
    reconstruction_count = options.nsamples
    sequence_length_string = options.sequence_length
    if sequence_length_string == 'inf':
        sequence_length = float('inf')
    else:
        sequence_length = int(sequence_length_string)
    inf_replacement = 20.0
    if options.reject_inf:
        inf_replacement = None
    elif options.replace_inf:
        try:
            inf_replacement = float(options.replace_inf)
        except ValueError:
            msg = 'invalid replace_inf value: '
            raise OptionError(msg + str(options.replace_inf))
    zero_replacement = 0
    if options.reject_zero:
        zero_replacement = None
    elif options.replace_zero:
        try:
            zero_replacement = float(options.replace_zero)
        except ValueError:
            msg = 'invalid replace_zero value: '
            raise OptionError(msg + str(options.replace_zero))
    # start the html file
    print '<html><body>'
    # show the simulation parameters
    print 'original tree source:', tree_remark, '<br/>'
    print 'reconstruction count:', reconstruction_count, '<br/>'
    print 'sequence length:', sequence_length, '<br/>'
    # set the simulation parameters for each simulation
    for sim in sims:
        sim.set_original_tree(tree)
        # If there is only one reconstruction per method
        # then show the progress of the tree builder.
        if reconstruction_count == 1:
            sim.set_verbose()
    # define an arbitrary but consistent ordering of the taxa
    ordered_names = [node.name for node in tree.gen_tips()]
    try:
        # attempt to simulate a bunch of distance matrices
        if options.verbose:
            print 'sampling', reconstruction_count, 'distance matrices...'
        # initialize the distance matrix sampler
        sampler = DistanceMatrixSampler.DistanceMatrixSampler(
                tree, ordered_names, sequence_length)
        sampler.set_inf_replacement(inf_replacement)
        sampler.set_zero_replacement(zero_replacement)
        # start the progress bar
        pbar = Progress.Bar(1.0)
        # sample some distance matrices
        distance_matrices = []
        for result in sampler.gen_samples_or_none():
            # if we got a result then update the distance matrix list
            if result:
                sequence_list, D = result
                distance_matrices.append(D)
            # Update the progressbar regardless of whether or not
            # the proposal was accepted.
            remaining_acceptances = reconstruction_count - len(
                    distance_matrices)
            numerator = sampler.get_completed_proposals()
            denominator = numerator + sampler.get_remaining_proposals(
                    remaining_acceptances)
            dms_fraction = float(numerator) / float(denominator)
            dms_total = 1.0 / (1 + len(sims))
            pbar.update(dms_fraction * dms_total)
            # if we have enough samples then break the loop
            if not remaining_acceptances:
                break
        # reconstruct trees using various methods
        for i, sim in enumerate(sims):
            if options.verbose:
                print 'running "%s"...' % sim.description
            sim.run(distance_matrices, ordered_names)
            pbar.update(float(i+2) / float(1 + len(sims)))
        # stop the progress bar
        pbar.finish()
        # get the simulation data
        table = [('method', 'seconds', 'uniform loss', 'weighted loss')]
        for sim in sims:
            table.append((sim.description, sim.get_running_time(),
                sim.get_uniform_loss(), sim.get_deep_loss()))
        # convert the row major matrix into an html table
        print HtmlTable.get_table_string(table)
        # end the html file
        print '</html></body>'
    except KeyboardInterrupt:
        print 'interrupted stage', pbar.progress, 'of', pbar.high

def main(options):
    """
    Run some tree reconstructions from the command line.
    """
    # assert that various options are compatible
    if options.reject_inf and options.replace_inf:
        msg = 'reject_inf and replace_inf are incompatible options'
        raise OptionError(msg)
    if options.reject_zero and options.replace_zero:
        msg = 'reject_zero and replace_zero are incompatible options'
        raise OptionError(msg)
    # See if we are supposed to do some hard coded analysis,
    # possibly using a user tree.
    if options.hard_coded_analysis:
        tree, tree_remark = get_tree_and_remark(options)
        if options.hard_coded_analysis == 1:
            do_hard_coded_analysis_a(tree, tree_remark)
        elif options.hard_coded_analysis == 2:
            do_hard_coded_analysis_b(tree, tree_remark)
        else:
            msg = 'invalid hard_coded_analysis: '
            raise OptionError(msg + options.hard_coded_analysis)
    else:
        do_command_line_analysis(options)

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('-v', '--verbose',
            action='store_true', dest='verbose',
            default=False, help='show extra information')
    parser.add_option('--nsamples', dest='nsamples', type='int', 
            default=1000, help='number of samples')
    parser.add_option('--sequence-length', dest='sequence_length',
            default=100,
            help='columns in the intermediate sampled alignment (or inf)')
    parser.add_option('--tree', dest='tree_filename',
            default='', help='tree filename')
    parser.add_option('--inline-tree', dest='inline_tree',
            default='', help='newick tree specified on the command line')
    parser.add_option('--use-exact', action='store_true', dest='use_exact',
            default=False, help='use the exact criterion (slow)')
    parser.add_option('--reject-inf', action='store_true', dest='reject_inf',
            default=False, help='reject matrices with distances of infinity')
    parser.add_option('--reject-zero', action='store_true', dest='reject_zero',
            default=False, help='reject matrices with distances of zero')
    parser.add_option('--replace-inf', dest='replace_inf',
            default='', help='use a given distance instead of infinity')
    parser.add_option('--replace-zero', dest='replace_zero',
            default='', help='use a given distance instead of zero')
    parser.add_option('--hard-coded-analysis-a',
            action="store_const", const=1, dest="hard_coded_analysis",
            default=None, help="make R files for analysis A")
    parser.add_option('--hard-coded-analysis-b',
            action="store_const", const=2, dest="hard_coded_analysis",
            default=None, help="make R files for analysis B")
    options, args = parser.parse_args()
    try:
        main(options)
    except OptionError, e:
        print e

