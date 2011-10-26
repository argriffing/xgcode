"""Create position-specific annotations of resequenced data.

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
A: an integer representing the number of A reads aligned to the this position
C: an integer representing the number of C reads aligned to the this position
G: an integer representing the number of G reads aligned to the this position
T: an integer representing the number of T reads aligned to the this position
hom_ll: the floating point non-Markov log likelihood
that the position is homozygous
het_ll: the floating point non-Markov log likelihood
that the position is heterozygous
bad_ll: the floating point non-Markov log likelihood
that the position is bad
hom_post: the floating point posterior probability
that the position is homozygous
het_post: the floating point posterior probability
that the position is heterozygous
bad_post: the floating point posterior probability
that the position is bad
"""

from StringIO import StringIO
import time
import sys

import argparse
import numpy as np

from SnippetUtil import HandlingError
import Form
import FormOut
import Progress
import ReadCoverage
import MissingHMM
import iterutils
import const

g_sample_data = const.read('20100730k')

g_output_header_row = [
        'genetic_line', 'chromosome', 'position',
        'A_count', 'C_count', 'G_count', 'T_count',
        'homozygous_ll', 'heterozygous_ll', 'bad_ll',
        'homozygous_posterior', 'heterozygous_posterior', 'bad_posterior']

class TimeoutError(Exception): pass

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.Integer('good_coverage',
                'expected read coverage of informative positions',
                20, low=1, high=1000),
            Form.Integer('bad_coverage',
                'expected read coverage of overcovered positions',
                100, low=1, high=1000),
            Form.Float('randomization_rate',
                'randomization probability per base call',
                0.1, low_exclusive=0),
            Form.MultiLine('input_text',
                'calls per nt per base call per chromosome per strain',
                g_sample_data),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.Csv()

def get_response_content(fs):
    # allow only two seconds for web access, and don't use a progress bar
    nseconds = 2
    use_pbar = False
    # get the lines from the multi-line input
    lines = StringIO(fs.input_text).readlines()
    lines = [line.strip() for line in lines]
    lines = [line for line in lines if line]
    # define an arbitrary transition matrix
    T = np.array([
        [.98, .01, .01],
        [.01, .98, .01],
        [.01, .01, .98]])
    # try to get the response
    try:
        response_text = process(lines, fs.good_coverage, fs.bad_coverage, fs.randomization_rate, T, nseconds, use_pbar)
    except TimeoutError:
        raise HandlingError('exceeded the remote time limit')
    # return the response
    return response_text + '\n'


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
        except ValueError as e
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
        @param coverages: the A, C, G, T coverage counts
        """
        if len(coverages) != 4:
            raise ValueError('four coverage counts were expected')
        self.offsets.append(offset)
        self.nt_coverages.append(coverages)

    def annotate_likelihoods(self, hidden_models):
        """
        @param hidden_models: a list of statistical models
        """
        observations = self.nt_coverages
        self.log_likelihoods = [[m.get_log_likelihood(obs) for m in hidden_models] for obs in observations]

    def annotate_posteriors(self, T, hidden_models):
        """
        @param T: a matrix of transition probabilities among the hidden states
        @param hidden_models: a list of statistical models
        """
        # define the HMM
        hmm = MissingHMM.MissingHMM(T, hidden_models)
        # define the observations and distances
        observations = self.nt_coverages
        distances = [b - a for a, b in iterutils.pairwise(self.offsets)]
        # do the annotation
        self.posterior_distributions = hmm.scaled_posterior_durbin(observations, distances)

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
    if len(input_lines) < 6:
        raise ValueError('there should be at least six lines of input')
    if len(input_lines) % 5 != 1:
        raise ValueError('the input lines should consist of a header plus a multiple of five data lines')
    # break the lines of input into rows of elements
    input_rows = [line_to_row(line) for line in input_lines]
    # validate the columns of data
    ncolumns = len(input_rows[0])
    if ncolumns < 7:
        raise ValueError('there should be at least seven columns of input')
    if ncolumns % 2 != 1:
        raise ValueError('the number of input columns should be odd')
    for row in input_rows:
        if len(row) != ncolumns:
            raise ValueError('each row of input should have the same number of elements as the first row')
    # read some data from the header
    header_row = input_rows[0]
    strain_names = []
    for heading in header_row[5:]:
        if heading.endswith('sco'):
            pass
        elif heading.endswith('cov'):
            strain_names.append(heading[:-3])
        else:
            raise ValueError('each heading after the fifth should end with sco or cov')
    # get the rest of the rows
    data_rows = input_rows[1:]
    # define the number of genomic positions and the number of strains
    npositions = len(data_rows) / 5
    nstrains = (ncolumns - 5) / 2
    # parse the data rows
    chromosome_dict = {}
    for position in range(npositions):
        # get a chunk of five consecutive rows
        position_rows = [data_rows[position*5 + i] for i in range(5)]
        # parse the chunk of five rows
        process_genomic_position(position_rows, strain_names, chromosome_dict)
    # turn the dictionary of chromosomes into a somewhat ordered list
    chromosomes = [chromosome for identifier, chromosome in sorted(chromosome_dict.items())]
    # return the list of chromosomes
    return chromosomes

def process_genomic_position(rows, strain_names, chromosome_dict):
    """
    The input data represents a single genomic position.
    Each chromosome object is accessible by its (strain_name, chromosome_name) pair.
    @param rows: five rows of the csv file
    @param strain_names: an ordered list of strain names
    @param chromosome_dict: a dictionary of chromosome objects
    """
    # Read the chromosome name from the first of the five rows.
    chromosome_name = rows[0][0]
    # Read the genomic position from the first of the five rows.
    genomic_position = int(rows[0][3])
    # Extract the important columns of the important rows,
    # and convert the elements from strings to integers.
    reduced_rows = []
    for row in rows[:4]:
        reduced_row = []
        for i, element in enumerate(row):
            if i > 4 and i % 2 == 1:
                count = int(element)
                if count < 0:
                    raise ValueError('counts must be nonnegative: ' + str(count))
                reduced_row.append(count)
        reduced_rows.append(reduced_row)
    # Each reduced row corresponds to a nucleotide.
    # Each column in a reduced row corresponds to a strain.
    # Each element is a coverage count.
    observations = zip(*reduced_rows)
    # Add an observation to each chromosome.
    for observation, strain_name in zip(observations, strain_names):
        # The input observation order is ACTG,
        # but the order I want to use is ACGT.
        A, C, T, G = observation
        observation = (A, C, G, T)
        # get or create the chromosome object
        identifier = (strain_name, chromosome_name)
        if identifier in chromosome_dict:
            chromosome = chromosome_dict[identifier]
        else:
            chromosome = Chromosome(*identifier)
            chromosome_dict[identifier] = chromosome
        # add the position to the chromosome
        chromosome.add_position(genomic_position, observation)

def line_to_row(line):
    """
    @param line: a line in the csv file
    @return: a list of comma separated elements
    """
    return [x.strip() for x in line.split(',')]

def process(input_lines, good_coverage, bad_coverage, randomization_rate, T, nseconds, use_pbar):
    """
    @param input_lines: lines of input of csv data including the header
    @param good_coverage: the expected number of reads at informative positions
    @param bad_coverage: the expected number of reads at uninformative positions
    @param randomization_rate: the probability of an error per read
    @param T: a transition matrix relating the hidden states
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
            chromosome.annotate_posteriors(T, models)
            print >> out, '\n'.join(','.join(row) for row in chromosome.get_rows_of_strings())
            if pbar:
                pbar.update(i + 1)
    except KeyboardInterrupt, e:
        if pbar:
            pbar.finish()
        raise e
    except TimeoutError as e
        if pbar:
            pbar.finish()
        raise e
    # return the output text
    return out.getvalue().strip()

def main(args):
    # validate the options
    assert 1 <= args.good_coverage
    assert 1 <= args.bad_coverage
    assert 0 < args.randomization_rate <= 1
    # define an arbitrary transition matrix
    T = np.array([
        [.98, .01, .01],
        [.01, .98, .01],
        [.01, .01, .98]])
    # read from standard input
    lines = sys.stdin.readlines()
    lines = [line.strip() for line in lines]
    lines = [line for line in lines if line]
    # show the result
    use_pbar = True
    nseconds = None
    print process(lines,
            args.good_coverage, args.bad_coverage,
            args.randomization_rate, T, nseconds, use_pbar)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--good-coverage', type='int',
            default=20, help='expected read coverage of informative positions')
    parser.add_argument('--bad-coverage', type='int',
            default=100, help='expected read coverage of overcovered positions')
    parser.add_argument('--randomization-rate', type='float',
            default=0.1, help='randomization probability per read')
    main(parser.parse_args())
