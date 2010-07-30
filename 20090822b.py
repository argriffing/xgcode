"""Get two-level posterior state expectations from resequencing data.
"""

from StringIO import StringIO
import time
import optparse
import sys
import profile
import math

import numpy as np

from SnippetUtil import HandlingError
import Form
import FormOut
import Progress
import ReadCoverageGap
import FastHMM
import TransitionMatrix
import iterutils
import Util
import const

g_sample_data = const.read('20100730j')

class TimeoutError(Exception): pass

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.Integer('good_coverage',
                'expected read coverage of informative positions',
                20, low=1, high=1000),
            Form.Float('randomization_rate',
                'randomization probability per base call',
                0.1, low_exclusive=0),
            Form.Integer('stickiness',
                'level of stickiness',
                4, low=1, high=4),
            Form.MultiLine('input_text',
                'calls per nt per base call per chromosome per strain',
                g_sample_data),
            Form.RadioGroup('delivery', 'delivery', [
                Form.RadioItem('inline', 'view as text', True),
                Form.RadioItem('attachment', 'download as an R table')])]
    return form_objects

def get_form_out():
    return FormOut.RTable()

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # allow only two seconds for web access, and don't use a progress bar
    nseconds = 2
    use_pbar = False
    # make the multiline input look like one of many open files
    linesources = [StringIO(fs.input_text)]
    # try to get the response
    try:
        response_text = process(linesources, fs.good_coverage, fs.randomization_rate, fs.stickiness, nseconds, use_pbar)
    except TimeoutError:
        raise HandlingError('sorry scripts run remotely have the attention span of a fruit fly')
    # deliver the response in the appropriate format
    response_headers = [('Content-Type', 'text/plain')]
    if fs.attachment:
        output_filename = 'posterior-expectations.table'
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
        # initialize lists that will be read from the input file
        self.offsets = []
        self.nt_coverages = []
        # initialize the counts
        self.expected_count_vectors = []
        # initialize some other stuff
        self.ntransitions_expected = None
        self.log_likelihood = None

    def get_identifier(self):
        return self.strain_name, self.chromosome_name

    def del_position_specific_data(self):
        del self.offsets
        del self.nt_coverages

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
        # use unlimited cache sizes
        cache_limit = None
        # define the transition matrix
        nhidden = len(hidden_models)
        prandom = .1**stickiness
        transition_object = TransitionMatrix.UniformTransitionObject(prandom, nhidden, cache_limit)
        # define the HMM
        hmm = FastHMM.Model(transition_object, hidden_models, cache_limit)
        # define the observations and distances
        observations = [tuple(sorted(coverage)) for coverage in self.nt_coverages]
        distances = [b - a for a, b in iterutils.pairwise(self.offsets)]
        # get the posterior distribution for each observation
        dp_info = hmm.get_dp_info(observations, distances)
        distribution_list = hmm.scaled_posterior_durbin(dp_info)
        # initialize the counts
        for model in hidden_models:
            self.expected_count_vectors.append(np.zeros(len(model.states)))
        # accumulate the counts
        for observation, distribution in zip(observations, distribution_list):
            for p in distribution:
                if math.isnan(p):
                    raise ValueError('nan in distribution: %s' % distribution)
            vectors = [model.get_posterior_distribution(observation) for model in hidden_models]
            for v in vectors:
                for x in v:
                    if math.isnan(x):
                        raise ValueError('nan in posterior mixture: %s' % v)
            normalized_vectors = [v*p for v, p in zip(vectors, distribution)]
            for i, v in enumerate(normalized_vectors):
                self.expected_count_vectors[i] += v
        # compute the log likelihood
        self.log_likelihood = hmm.get_log_likelihood(dp_info)
        # compute the expected number of hidden state transitions
        self.ntransitions_expected = hmm.scaled_ntransitions_expected(dp_info)

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

def process(linesources, good_coverage, randomization_rate, stickiness, nseconds, use_pbar):
    """
    @param linesources: open resequencing files for reading
    @param good_coverage: the expected number of reads at informative positions
    @param randomization_rate: the probability of an error per base call
    @param stickiness: level of stickiness
    @param nseconds: None or impose a time limit of this many seconds
    @param use_pbar: True iff a progress bar should be used
    @return: the multi-line string of the resulting csv file
    """
    # do some initialization
    start_time = time.time()
    termination_reason = 'finished the analysis'
    # define the superstates
    cache_size_per_superstate = 100000
    good_state = ReadCoverageGap.Good(randomization_rate, good_coverage, cache_size_per_superstate)
    bad_state = ReadCoverageGap.Bad(randomization_rate, good_coverage, cache_size_per_superstate)
    superstates = [good_state, bad_state]
    superstate_names = ['good', 'bad']
    # prepare to annotate the chromosomes
    chromosomes = []
    pbar = Progress.Bar(len(linesources)) if use_pbar else None
    # annotate the chromosomes using the models
    try:
        for i, linesource in enumerate(linesources):
            # read the lines of text
            lines = Util.get_stripped_lines(linesource.readlines())
            # validate the number of lines
            if len(lines) < 2:
                raise ValueError('there should be at least two lines of input')
            # break the lines of input into rows of elements
            rows = [line_to_row(line) for line in lines]
            # validate the columns of data
            ncolumns_expected = 8
            ncolumns = len(rows[0])
            if ncolumns != ncolumns_expected:
                raise ValueError('expected %d columns of input: %s' % (ncolumns_expected, rows[0]))
            for row in rows:
                if len(row) != ncolumns:
                    raise ValueError('each row of input should have the same number of elements as the first row')
            # process the data rows, building a dictionary of chromosomes
            chromosome_dict = {}
            data_rows = rows[1:]
            for row in data_rows:
                if nseconds and time.time() - start_time > nseconds:
                    raise TimeoutError()
                process_genomic_position(row, chromosome_dict)
            current_chromosomes = [chromosome for identifier, chromosome in sorted(chromosome_dict.items())]
            for chromosome in current_chromosomes:
                # do the annotation
                chromosome.annotate_posteriors(stickiness, superstates)
                # delete position specific data
                chromosome.del_position_specific_data()
            # add the chromosomes to the list
            chromosomes.extend(current_chromosomes)
            # update the progress bar
            if pbar:
                pbar.update(i + 1)
    except KeyboardInterrupt, e:
        termination_reason = 'early termination by control-c'
    except TimeoutError, e:
        termination_reason = 'early termination because of a time limit'
    if pbar:
        pbar.finish()
    # begin the response
    out = StringIO()
    # write some meta data
    print >> out, '#', 'termination:', termination_reason
    print >> out, '#', 'elapsed seconds:', time.time() - start_time
    print >> out, '#', 'good coverage:', good_coverage
    print >> out, '#', 'stickiness:', stickiness
    print >> out, '#', 'randomization_rate:', randomization_rate
    # write the header line
    header_row = ['strain', 'chromosome']
    for superstate, superstate_name in zip(superstates, superstate_names):
        for substate_name in superstate.state_names:
            heading = '%s.%s' % (superstate_name, substate_name)
            header_row.append(heading)
    header_row.extend(['log.likelihood', 'expected.transitions'])
    print >> out, '\t'.join(header_row)
    # reorder the list of chromosomes
    identifier_chromosome_pairs = [(x.get_identifier(), x) for x in chromosomes]
    chromosomes = [chrom for identifier, chrom in sorted(identifier_chromosome_pairs)]
    # write the data lines
    for i, chromosome in enumerate(chromosomes):
        row = []
        row.extend([str(i), chromosome.strain_name, chromosome.chromosome_name])
        for v in chromosome.expected_count_vectors:
            row.extend([str(x) for x in v.tolist()])
        row.extend([str(chromosome.log_likelihood), str(chromosome.ntransitions_expected)])
        print >> out, '\t'.join(row)
    # return the result
    return out.getvalue().strip()

def main(options, pathnames):
    # validate the options
    assert 1 <= options.good_coverage
    assert 0 < options.randomization_rate <= 1
    assert 0 <= options.stickiness
    # show the result
    use_pbar = True
    nseconds = None
    linesources = [open(pathname) for pathname in pathnames]
    print process(linesources, options.good_coverage, options.randomization_rate, options.stickiness, nseconds, use_pbar)
    for linesource in linesources:
        linesource.close()

if __name__ == '__main__':
    from optparse import OptionParser
    usage = 'Usage: %prog <file_1.csv> <file_2.csv> ... <file_n.csv>'
    parser = OptionParser(usage=usage)
    parser.add_option('--good-coverage', dest='good_coverage', type='int', default=20, help='expected read coverage of informative positions')
    parser.add_option('--randomization-rate', dest='randomization_rate', type='float', default=0.1, help='randomization probability per base call')
    parser.add_option('--stickiness', dest='stickiness', type='int', default=4, help='level of stickiness')
    parser.add_option('--profile', action='store_true', dest='profile')
    options, args = parser.parse_args()
    if options.profile:
        profile.run('main(options, args)')
    else:
        main(options, args)
