"""Get some more statistics for resequenced chromosomes.

Note that this analysis does *not* use a Markov model.
For each input file representing a chromosome of a genetic line,
compute the expected amount of time spent in each of three states.
Columns of each non-header row of each input file should be as follows.
strain: a string defining the genetic line
chromosome: a string defining the name of the chromosome
offset: an integer offset in units of base pairs
A: an integer representing the number of A reads aligned to the this position
C: an integer representing the number of C reads aligned to the this position
G: an integer representing the number of G reads aligned to the this position
T: an integer representing the number of T reads aligned to the this position
gap: an integer representing the number of gaps reads
aligned to the this position
"""

from StringIO import StringIO
import time
import optparse
import sys
import math

import numpy as np

from SnippetUtil import HandlingError
import Form
import Progress
import ReadCoverage


class TimeoutError(Exception): pass


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


class Oracle:
    def __init__(self, models, cache_size):
        """
        @param models: get a likelihood from a model
        @param cache_size: the number of observations to remember
        """
        self.models = models
        self.cache_size = cache_size
        self.cache = {}
    def get_model_count(self):
        return len(self.models)
    def get_distribution(self, obs):
        """
        This function uses memoization.
        @param obs: an emitted state
        @return: a stochastic vector as a numpy array
        """
        v = self.cache.get(obs, None)
        if v is not None:
            return v
        v = np.array([m.get_log_likelihood(obs) for m in self.models])
        v = np.exp(v - np.max(v))
        v /= np.sum(v)
        if len(self.cache) < self.cache_size:
            self.cache[obs] = v
        return v


class Chromosome:
    """
    One of these exists for each chromosome in each strain.
    """

    def __init__(self, strain_name, chromosome_name, oracle):
        """
        @param strain_name: the name of the associated genetic strain
        @param chromosome_name: the name of the chromosome
        @param oracle: a source of distribution vectors
        """
        try:
            test_string = ','.join([strain_name, chromosome_name])
        except ValueError, e:
            raise ValueError('value error for names %s and %s: %s' % (strain_name, chromosome_name, e))
        # initialize some states
        self.strain_name = strain_name
        self.chromosome_name = chromosome_name
        self.oracle = oracle
        # initialize stuff
        self.expected_counts = np.zeros(self.oracle.get_model_count())

    def add_position(self, coverage):
        """
        Add a position to the chromosome.
        @param coverage: the (A, C, G, T, gap) coverage counts
        """
        if len(coverage) != 5:
            raise ValueError('five coverage counts were expected')
        obs = tuple(sorted(coverage[:-1]))
        self.expected_counts += self.oracle.get_distribution(obs)


def process_genomic_position(row, chromosome_dict, oracle):
    """
    The input data represents a single genomic position.
    Each chromosome object is accessible by its (strain_name, chromosome_name) pair.
    @param row: a data row from the csv file
    @param chromosome_dict: a dictionary of chromosome objects
    @param oracle: has models
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
        chromosome = Chromosome(strain_name, chromosome_name, oracle)
        chromosome_dict[identifier] = chromosome
    # add the observation to the chromosome
    chromosome.add_position(coverages)

def line_to_row(line):
    """
    @param line: a line in the csv file
    @return: a list of comma separated elements
    """
    return [x.strip() for x in line.split(',')]

def process(pathnames, low_coverage, good_coverage, bad_coverage, randomization_rate, nseconds, use_pbar):
    """
    @param pathnames: paths to files to process
    @param low_coverage: the expected number of reads at undercovered positions
    @param good_coverage: the expected number of reads at informative positions
    @param bad_coverage: the expected number of reads at overcovered positions
    @param randomization_rate: the probability of an error per read
    @param nseconds: None or impose a time limit of this many seconds
    @param use_pbar: True iff a progress bar should be used
    @return: the multi-line string of the resulting csv file
    """
    # define the models
    homozygous = ReadCoverage.Homozygous(randomization_rate, good_coverage)
    heterozygous = ReadCoverage.Heterozygous(randomization_rate, good_coverage)
    overcovered = ReadCoverage.Overcovered(randomization_rate, bad_coverage)
    undercovered = ReadCoverage.FlatState(low_coverage)
    models = [homozygous, heterozygous, undercovered, overcovered]
    # define the oracle
    cache_size = 100000
    oracle = Oracle(models, cache_size)
    # do some initialization
    out = StringIO()
    start_time = time.time()
    nfiles = len(pathnames)
    pbar = Progress.Bar(nfiles) if (use_pbar and nfiles > 1) else None
    chromosome_dict = {}
    termination_reason = 'finished the analysis'
    try:
        for i, pathname in enumerate(pathnames):
            with open(pathname) as fin:
                lines = fin.readlines()
                lines = [line.strip() for line in lines]
                lines = [line for line in lines if line]
                # validate the number of lines
                if len(lines) < 2:
                    raise ValueError('there should be at least two lines of input')
                # break the lines of input into rows of elements
                rows = [line_to_row(line) for line in lines]
                # validate the columns of data
                ncolumns_expected = 8
                ncolumns = len(rows[0])
                if ncolumns != ncolumns_expected:
                    raise ValueError('expected %d columns of input' % ncolumns_expected)
                for row in rows:
                    if len(row) != ncolumns:
                        raise ValueError('each row of input should have the same number of elements as the first row')
                # process the data rows
                data_rows = rows[1:]
                for row in data_rows:
                    if nseconds and time.time() - start_time > nseconds:
                        raise TimeoutError()
                    process_genomic_position(row, chromosome_dict, oracle)
            if pbar:
                pbar.update(i + 1)
    except KeyboardInterrupt, e:
        termination_reason = 'early termination by control-c'
    except TimeoutError, e:
        termination_reason = 'early termination because of a time limit'
    if pbar:
        pbar.finish()
    # write some meta data
    print >> out, '#', 'note: this analysis does not use a markov model'
    print >> out, '#', 'termination:', termination_reason
    print >> out, '#', 'sys.argv:', ' '.join(sys.argv)
    print >> out, '#', 'elapsed seconds:', time.time() - start_time
    print >> out, '#', 'low_coverage:', low_coverage
    print >> out, '#', 'good_coverage:', good_coverage
    print >> out, '#', 'bad_coverage:', bad_coverage
    print >> out, '#', 'randomization_rate:', randomization_rate
    # write the header line
    print >> out, '\t'.join(['strain', 'chromosome', 'hom', 'het', 'low', 'other'])
    # turn the dictionary of chromosomes into a somewhat ordered list
    chromosomes = [chromosome for identifier, chromosome in sorted(chromosome_dict.items())]
    # write the data lines
    for i, chromosome in enumerate(chromosomes):
        row = []
        row.extend([str(i), chromosome.strain_name, chromosome.chromosome_name])
        row.extend([str(c_hat) for c_hat in chromosome.expected_counts])
        print >> out, '\t'.join(row)
    # return the result
    return out.getvalue().strip()

def main(options, pathnames):
    # validate the options
    assert 1 <= options.low_coverage
    assert 1 <= options.good_coverage
    assert 1 <= options.bad_coverage
    assert 0 < options.randomization_rate <= 1
    # show the result
    use_pbar = True
    nseconds = None
    print process(pathnames, options.low_coverage, options.good_coverage, options.bad_coverage, options.randomization_rate, nseconds, use_pbar)

if __name__ == '__main__':
    from optparse import OptionParser
    usage = 'Usage: %prog <file_1.csv> <file_2.csv> ... <file_n.csv>'
    parser = OptionParser(usage=usage)
    parser.add_option('--low-coverage', dest='low_coverage', type='int', default=2, help='expected read coverage of undercovered positions')
    parser.add_option('--good-coverage', dest='good_coverage', type='int', default=20, help='expected read coverage of informative positions')
    parser.add_option('--bad-coverage', dest='bad_coverage', type='int', default=100, help='expected read coverage of overcovered positions')
    parser.add_option('--randomization-rate', dest='randomization_rate', type='float', default=0.1, help='randomization probability per read')
    options, args = parser.parse_args()
    main(options, args)

