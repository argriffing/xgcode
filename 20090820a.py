"""Annotate a single resequenced chromosome in a single genetic line.

Use a HMM with 3 hidden states.
A subset of genomic positions are annotated with a posterior hidden state distribution
given observations at each position in the subset.
The 3x3 transition matrix is known,
and the emission distributions for each hidden state are known
but are somewhat complicated.
Columns of the output are defined as follows.
strain: a string defining the genetic line
chromosome: a string defining the name of the chromosome
offset: an integer offset in units of base pairs
A: an integer representing the number of A reads aligned to the this position
C: an integer representing the number of C reads aligned to the this position
G: an integer representing the number of G reads aligned to the this position
T: an integer representing the number of T reads aligned to the this position
gap: an integer representing the number of gaps reads aligned to the this position
hom_ll: the floating point non-Markov log likelihood that the position is homozygous
het_ll: the floating point non-Markov log likelihood that the position is heterozygous
bad_ll: the floating point non-Markov log likelihood that the position is bad
hom_post: the floating point posterior probability that the position is homozygous
het_post: the floating point posterior probability that the position is heterozygous
bad_post: the floating point posterior probability that the position is bad
"""

import StringIO
import time
import optparse
import sys
import profile

import numpy as np

from SnippetUtil import HandlingError
import Form
import Progress
import ReadCoverage
import FastHMM
import TransitionMatrix
import Util


class TimeoutError(Exception): pass


g_output_header_row = [
        'genetic_line', 'chromosome', 'position',
        'A_count', 'C_count', 'G_count', 'T_count', 'gap_count',
        'homozygous_ll', 'heterozygous_ll', 'bad_ll',
        'homozygous_posterior', 'heterozygous_posterior', 'bad_posterior']


def get_form():
    """
    @return: the body of a form
    """
    return []

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # allow only two seconds for web access, and don't use a progress bar
    nseconds = 2
    use_pbar = False
    # try to get the response
    try:
        response_text = 'for now this script is accessible only from the command line'
    except TimeoutError:
        response_text = 'sorry scripts run remotely have the attention span of a fruit fly'
    return [('Content-Type', 'text/plain')], response_text


class Chromosome:
    """
    One of these exists for each chromosome in each strain.
    This is the level at which the HMM applies.
    """

    def __init__(self, strain_name, chromosome_name):
        """
        @param strain_name: the name of the associated genetic strain
        @param chromosome_name: the name of the chromosome
        """
        try:
            test_string = ','.join([strain_name, chromosome_name])
        except ValueError, e:
            raise ValueError('value error for names %s and %s: %s' % (strain_name, chromosome_name, e))
        # initialize some states
        self.strain_name = strain_name
        self.chromosome_name = chromosome_name
        # initialize stuff
        self.offsets = []
        self.nt_coverages = []
        self.log_likelihoods = []
        self.posterior_distributions = []

    def add_position(self, offset, coverages):
        """
        Add a position to the chromosome.
        @param offset: the nucleotide offset
        @param coverages: the (A, C, G, T, gap) coverage counts
        """
        if len(coverages) != 5:
            raise ValueError('five coverage counts were expected')
        self.offsets.append(offset)
        self.nt_coverages.append(coverages)

    def annotate_likelihoods(self, hidden_models):
        """
        @param hidden_models: a list of statistical models
        """
        observations = [coverage[:-1] for coverage in self.nt_coverages]
        self.log_likelihoods = [[m.get_log_likelihood(obs) for m in hidden_models] for obs in observations]

    def annotate_posteriors(self, transition_object, hidden_models):
        """
        @param transition_object: has transition matrix information
        @param hidden_models: a list of statistical models
        """
        # define the HMM
        cache_size = 10000
        hmm = FastHMM.Model(transition_object, hidden_models, cache_size)
        # define the observations and distances
        observations = [tuple(sorted(coverage[:-1])) for coverage in self.nt_coverages]
        distances = [b - a for a, b in Util.pairwise(self.offsets)]
        # do the annotation
        dp_info = hmm.get_dp_info(observations, distances)
        self.posterior_distributions = hmm.scaled_posterior_durbin(dp_info)

    def get_rows_of_strings(self):
        """
        @return: a list of lists of strings
        """
        rows = []
        for i, offset in enumerate(self.offsets):
            row = [self.strain_name, self.chromosome_name, str(offset)]
            row.extend([str(x) for x in self.nt_coverages[i]])
            row.extend([str(x) for x in self.log_likelihoods[i]])
            row.extend([str(x) for x in self.posterior_distributions[i]])
            rows.append(row)
        return rows


def parse(input_lines):
    """
    @param input_lines: raw input lines of the csv file
    @return: a list of chromosome objects
    """
    # validate the number of lines
    if len(input_lines) < 2:
        raise ValueError('there should be at least two lines of input')
    # break the lines of input into rows of elements
    input_rows = [line_to_row(line) for line in input_lines]
    # validate the columns of data
    ncolumns_expected = 8
    ncolumns = len(input_rows[0])
    if ncolumns != ncolumns_expected:
        raise ValueError('expected %d columns of input' % ncolumns_expected)
    for row in input_rows:
        if len(row) != ncolumns:
            raise ValueError('each row of input should have the same number of elements as the first row')
    # process the non-header rows
    chromosome_dict = {}
    for row in input_rows[1:]:
        process_genomic_position(row, chromosome_dict)
    # turn the dictionary of chromosomes into a somewhat ordered list
    chromosomes = [chromosome for identifier, chromosome in sorted(chromosome_dict.items())]
    # return the list of chromosomes
    return chromosomes

def process_genomic_position(row, chromosome_dict):
    """
    The input data represents a single genomic position.
    Each chromosome object is accessible by its (strain_name, chromosome_name) pair.
    @param row: a data row from the csv file
    @param chromosome_dict: a dictionary of chromosome objects
    """
    # unpack the row
    strain_name, chromosome_name = row[:2]
    genomic_position = int(row[2])
    coverages = [int(x) for x in row[3:]]
    # get or create the relevant chromosome object
    identifier = (strain_name, chromosome_name)
    if identifier in chromosome_dict:
        chromosome = chromosome_dict[identifier]
    else:
        chromosome = Chromosome(*identifier)
        chromosome_dict[identifier] = chromosome
    # add the observation to the chromosome
    chromosome.add_position(genomic_position, coverages)

def line_to_row(line):
    """
    @param line: a line in the csv file
    @return: a list of comma separated elements
    """
    return [x.strip() for x in line.split(',')]

def process(input_lines, good_coverage, bad_coverage, randomization_rate, transition_object, nseconds, use_pbar):
    """
    @param input_lines: lines of input of csv data including the header
    @param good_coverage: the expected number of reads at informative positions
    @param bad_coverage: the expected number of reads at uninformative positions
    @param randomization_rate: the probability of an error per read
    @param transition_object: has transition matrix information
    @param nseconds: None or impose a time limit of this many seconds
    @param use_pbar: True iff a progress bar should be used
    @return: the multi-line string of the resulting csv file
    """
    # do some initialization
    out = StringIO.StringIO()
    pbar = None
    start_time = time.time()
    # define the three models
    homozygous = ReadCoverage.Homozygous(randomization_rate, good_coverage)
    heterozygous = ReadCoverage.Heterozygous(randomization_rate, good_coverage)
    overcovered = ReadCoverage.Overcovered(randomization_rate, bad_coverage)
    models = [homozygous, heterozygous, overcovered]
    # read the chromosome data
    chromosomes = parse(input_lines)
    # write the header line
    print >> out, ','.join(g_output_header_row)
    # prepare to annotate the chromosomes
    if use_pbar:
        pbar = Progress.Bar(len(chromosomes))
    # annotate the chromosomes using the models
    try:
        for i, chromosome in enumerate(chromosomes):
            if nseconds and time.time() - start_time > nseconds:
                raise TimeoutError()
            chromosome.annotate_likelihoods(models)
            chromosome.annotate_posteriors(transition_object, models)
            print >> out, '\n'.join(','.join(row) for row in chromosome.get_rows_of_strings())
            if pbar:
                pbar.update(i + 1)
    except KeyboardInterrupt, e:
        if pbar:
            pbar.finish()
        raise e
    except TimeoutError, e:
        if pbar:
            pbar.finish()
        raise e
    # return the output text
    return out.getvalue().strip()

def main(options):
    # validate the options
    assert 1 <= options.good_coverage
    assert 1 <= options.bad_coverage
    assert 0 < options.randomization_rate <= 1
    # define an arbitrary transition matrix
    transition_object = TransitionMatrix.UniformTransitionObject(.01, 3)
    # read from standard input
    lines = sys.stdin.readlines()
    lines = [line.strip() for line in lines]
    lines = [line for line in lines if line]
    # show the result
    use_pbar = True
    nseconds = None
    print process(lines, options.good_coverage, options.bad_coverage, options.randomization_rate, transition_object, nseconds, use_pbar)

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('--good-coverage', dest='good_coverage', type='int', default=20, help='expected read coverage of informative positions')
    parser.add_option('--bad-coverage', dest='bad_coverage', type='int', default=100, help='expected read coverage of overcovered positions')
    parser.add_option('--randomization-rate', dest='randomization_rate', type='float', default=0.1, help='randomization probability per read')
    parser.add_option('--profile', action='store_true', dest='profile')
    options, args = parser.parse_args()
    if options.profile:
        profile.run('main(options)')
    else:
        main(options)

