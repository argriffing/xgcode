"""Annotate a chromosome with MAP states for various prior autocorrelations.

This should use HMMs with various levels of 'stickiness' or 'autocorrelation'.
Technically, the transition matrix at each state has probability (1-2x/3) of
no transition in one step, and a probability of (x/3)
of each of the two possible transitions,
where x=(1/10)**alpha where alpha is the 'stickiness'
which is a nonnegative integer.
A stickiness of zero means that the position specific posterior distributions
do not reflect any prior belief of sequential dependence among hidden states.
Use a HMM with 3 hidden states.
A subset of genomic positions are annotated
with a posterior hidden state distribution
given observations at each position in the subset.
The 3x3 transition matrix is known,
and the emission distributions for each hidden state are known
but are somewhat complicated.
Columns of the output are defined as follows.
strain: a string defining the genetic line
chromosome: a string defining the name of the chromosome
offset: an integer offset in units of base pairs
A: the number of A reads aligned to the this position
C: the number of C reads aligned to the this position
G: the number of G reads aligned to the this position
T: the number of T reads aligned to the this position
gap: the number of gaps reads aligned to the this position
call_0: MAP state for no prior stickiness
hom_0: relative posterior belief of homozygosity
het_0: relative posterior belief of heterozygosity
bad_0: relative posterior belief of badness
call_1: MAP state for prior stickiness of 1
hom_1: relative posterior belief of homozygosity
het_1: relative posterior belief of heterozygosity
bad_1: relative posterior belief of badness
.
.
.
"""

from StringIO import StringIO
import time
import optparse
import sys
import profile

import numpy as np

from SnippetUtil import HandlingError
import Form
import Progress
import TransitionMatrix
import ReadCoverage
import FastHMM
import iterutils


class TimeoutError(Exception): pass


g_sample_rows = [
        ['strain', 'chromosome', 'position', 'A', 'C', 'G', 'T', 'gap'],
        ['222', '1024', '101', '0', '0', '0', '0', '0'],
        ['222', '1024', '121', '21', '1', '0', '0', '0'],
        ['222', '1024', '151', '0', '12', '2', '0', '9'],
        ['444', '1024', '201', '0', '0', '0', '0', '0'],
        ['444', '1024', '231', '21', '1', '0', '0', '0'],
        ['444', '1024', '241', '0', '12', '2', '666', '9']]

def get_form():
    """
    @return: the body of a form
    """
    sample_lines = [',\t'.join(row) for row in g_sample_rows]
    form_objects = [
            Form.Integer('good_coverage',
                'expected read coverage of informative positions',
                20, low=1, high=1000),
            Form.Integer('bad_coverage',
                'expected read coverage of overcovered positions',
                100, low=1, high=1000),
            Form.Integer('nstickinesses',
                'use this many different levels of stickiness',
                5, low=1, high=5),
            Form.Float('randomization_rate',
                'randomization probability per base call',
                0.1, low_exclusive=0),
            Form.MultiLine('input_text',
                'calls per nt per base call per chromosome per strain',
                '\n'.join(sample_lines)),
            Form.RadioGroup('delivery', 'delivery', [
                Form.RadioItem('inline', 'view as text', True),
                Form.RadioItem('attachment', 'download as a csv file')])]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # allow only two seconds for web access, and don't use a progress bar
    nseconds = 2
    use_pbar = False
    # get the lines from the multi-line input
    lines = StringIO(fs.input_text).readlines()
    lines = [line.strip() for line in lines]
    lines = [line for line in lines if line]
    # try to get the response
    try:
        response_text = process(lines, fs.good_coverage, fs.bad_coverage, fs.randomization_rate, fs.nstickinesses, nseconds, use_pbar)
    except TimeoutError:
        raise HandlingError('sorry scripts run remotely have the attention span of a fruit fly')
    # deliver the response in the appropriate format
    response_headers = [('Content-Type', 'text/plain')]
    if fs.attachment:
        output_filename = 'annotation.csv'
        response_headers.append(('Content-Disposition', "%s; filename=%s" % (fs.delivery, output_filename)))
    return response_headers, response_text


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
        self.posterior_distribution_lists = []
        self.stickinesses = []

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

    def annotate_posteriors(self, stickiness, hidden_models):
        """
        @param stickiness: a nonnegative integer that defines the transition matrix
        @param hidden_models: a list of statistical models
        """
        # define the transition matrix
        nhidden = len(hidden_models)
        prandom = .1**stickiness
        transition_object = TransitionMatrix.UniformTransitionObject(prandom, nhidden)
        # define the HMM
        cache_size = 10000
        hmm = FastHMM.Model(transition_object, hidden_models, cache_size)
        # define the observations and distances
        observations = [tuple(sorted(coverage[:-1])) for coverage in self.nt_coverages]
        distances = [b - a for a, b in iterutils.pairwise(self.offsets)]
        # do the annotation
        dp_info = hmm.get_dp_info(observations, distances)
        distribution_list = hmm.scaled_posterior_durbin(dp_info)
        # store the annotation with its respective stickiness
        self.posterior_distribution_lists.append(distribution_list)
        self.stickinesses.append(stickiness)

    def get_rows_of_strings(self, model_names):
        """
        @param model_names: a conformant list of strings for annotation
        @return: a list of lists of strings
        """
        rows = []
        for i, offset in enumerate(self.offsets):
            row = []
            # append a fixed number of columns of non-annotation
            row.extend([self.strain_name, self.chromosome_name, str(offset)])
            row.extend([str(x) for x in self.nt_coverages[i]])
            # append a variable number columns of annotation
            for stickinesses, posterior_distributions in zip(self.stickinesses, self.posterior_distribution_lists):
                # identify the position specific posterior distribution for this stickiness
                distribution = posterior_distributions[i]
                # append the name of the MAP state defined by the distribution
                p, name = max(zip(distribution, model_names))
                row.append(name)
                # append the distribution itself
                row.extend([str(x) for x in distribution])
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

def process(input_lines, good_coverage, bad_coverage, randomization_rate, nstickinesses, nseconds, use_pbar):
    """
    @param input_lines: lines of input of csv data including the header
    @param good_coverage: the expected number of reads at informative positions
    @param bad_coverage: the expected number of reads at uninformative positions
    @param randomization_rate: the probability of an error per read
    @param nstickinesses: use this many different levels of stickiness
    @param nseconds: None or impose a time limit of this many seconds
    @param use_pbar: True iff a progress bar should be used
    @return: the multi-line string of the resulting csv file
    """
    # do some initialization
    out = StringIO()
    pbar = None
    start_time = time.time()
    # define the three models
    homozygous = ReadCoverage.Homozygous(randomization_rate, good_coverage)
    heterozygous = ReadCoverage.Heterozygous(randomization_rate, good_coverage)
    overcovered = ReadCoverage.Overcovered(randomization_rate, bad_coverage)
    models = [homozygous, heterozygous, overcovered]
    model_names = ['hom', 'het', 'ovr']
    # read the chromosome data
    chromosomes = parse(input_lines)
    # write the header line
    header_row = []
    header_row.extend([
        'genetic_line', 'chromosome', 'position',
        'A_count', 'C_count', 'G_count', 'T_count', 'gap_count'])
    for stickiness in range(nstickinesses):
        for name in ('call', 'hom', 'het', 'bad'):
            header_row.append('%s_%d' % (name, stickiness))
    print >> out, ','.join(header_row)
    # prepare to annotate the chromosomes
    if use_pbar:
        count = 0
        pbar = Progress.Bar(len(chromosomes)*nstickinesses)
    # annotate the chromosomes using the models
    for i, chromosome in enumerate(chromosomes):
        for stickiness in range(nstickinesses):
            if nseconds and time.time() - start_time > nseconds:
                raise TimeoutError()
            chromosome.annotate_posteriors(stickiness, models)
            if pbar:
                count += 1
                pbar.update(count)
        print >> out, '\n'.join(','.join(row) for row in chromosome.get_rows_of_strings(model_names))
    if pbar:
        pbar.finish()
    # return the output text
    return out.getvalue().strip()

def main(options):
    # validate the options
    assert 1 <= options.good_coverage
    assert 1 <= options.bad_coverage
    assert 0 < options.randomization_rate <= 1
    # read from standard input
    lines = sys.stdin.readlines()
    lines = [line.strip() for line in lines]
    lines = [line for line in lines if line]
    # show the result
    use_pbar = True
    nseconds = None
    print process(lines, options.good_coverage, options.bad_coverage, options.randomization_rate, options.nstickinesses, nseconds, use_pbar)

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('--nstickinesses', dest='nstickinesses', type='int', default=5, help='use this many different levels of stickiness')
    parser.add_option('--good-coverage', dest='good_coverage', type='int', default=20, help='expected read coverage of informative positions')
    parser.add_option('--bad-coverage', dest='bad_coverage', type='int', default=100, help='expected read coverage of overcovered positions')
    parser.add_option('--randomization-rate', dest='randomization_rate', type='float', default=0.1, help='randomization probability per read')
    parser.add_option('--profile', action='store_true', dest='profile')
    options, args = parser.parse_args()
    if options.profile:
        profile.run('main(options)')
    else:
        main(options)

