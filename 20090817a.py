"""Add log likelihood annotation to a csv file of genomic read counts.
"""

# Here is some terminology used in this file.
# strain: a genetic line of Drosophila, for example
# position: a genomic position including a chromosome and nucleotide offset
# line: a single line of unparsed data from a file
# row: a list of elements parsed from a comma separated line

from StringIO import StringIO
import time
import optparse
import sys

from SnippetUtil import HandlingError
import Form
import Progress
import ReadCoverage


class TimeoutError(Exception): pass


g_example_lines = [
        'chr,ref,var,pos,bp,208cov,208sco,301cov,301sco',
        'dmel_mitochondrion_genome,C,T,671,A,150,0,610,0',
        'dmel_mitochondrion_genome,C,T,671,C,12,2,2,314',
        'dmel_mitochondrion_genome,C,T,671,T,2,4,0,1337',
        'dmel_mitochondrion_genome,C,T,671,G,1,6,0,1',
        'dmel_mitochondrion_genome,C,T,671,-,0,0,0,0',
        'dmel_mitochondrion_genome,G,A,710,A,10,10,0,2',
        'dmel_mitochondrion_genome,G,A,710,C,7,12,0,3',
        'dmel_mitochondrion_genome,G,A,710,T,0,14,0,9',
        'dmel_mitochondrion_genome,G,A,710,G,0,16,0,81',
        'dmel_mitochondrion_genome,G,A,710,-,0,0,0,0',
        'dmel_mitochondrion_genome,G,A,710,A,2,10,19,2',
        'dmel_mitochondrion_genome,G,A,710,C,0,12,1,3',
        'dmel_mitochondrion_genome,G,A,710,T,0,14,0,9',
        'dmel_mitochondrion_genome,G,A,710,G,0,16,0,81',
        'dmel_mitochondrion_genome,G,A,710,-,0,0,0,0']

def get_form():
    """
    @return: the body of a form
    """
    # define the objects
    form_objects = [
            Form.Integer('good_coverage', 'expected coverage of homozygous and heterozygous positions', 20, low=1),
            Form.Integer('bad_coverage', 'expected coverage of overcovered positions', 100, low=1),
            Form.Float('randomization_rate', 'read randomization rate', 0.1, low_exclusive=0, high_inclusive=1),
            Form.MultiLine('csv_data', 'comma separated values', '\n'.join(g_example_lines))]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # allow only two seconds for web access, and don't use a progress bar
    nseconds = 2
    use_pbar = False
    # read the simple options
    good_coverage = fs.good_coverage
    bad_coverage = fs.bad_coverage
    randomization_rate = fs.randomization_rate
    # read the lines of comma separated values
    lines = StringIO(fs.csv_data).readlines()
    lines = [line.strip() for line in lines]
    lines = [line for line in lines if line]
    # try to get the response
    try:
        response_text = process(lines, good_coverage, bad_coverage, randomization_rate, nseconds, use_pbar)
    except TimeoutError:
        response_text = 'sorry scripts run remotely have the attention span of a fruit fly'
    return [('Content-Type', 'text/plain')], response_text

def get_log_likelihoods_per_strain(rows, models):
    """
    Return a list of log likelihood lists for each strain.
    Each log likelihood list is a list of log likelihoods of three models.
    The input data represents a single genomic position.
    @param rows: five rows of the csv file
    @param models: three models that can give a log likelihood from an observation
    @return: a list of likelihood lists
    """
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
    # get the log likelihoods
    log_likelihood_lists = []
    for obs in observations:
        assert len(obs) == 4
        n = sum(obs)
        if n > 0:
            try:
                log_likelihoods = [m.get_log_likelihood(obs) for m in models]
            except ValueError, e:
                raise ValueError('error for observation %s: %s' % (str(obs), e))
        else:
            log_likelihoods = None
        log_likelihood_lists.append(log_likelihoods)
    return log_likelihood_lists

def line_to_row(line):
    """
    @param line: a line in the csv file
    @return: a list of comma separated elements
    """
    return [x.strip() for x in line.split(',')]

def process(input_lines, good_coverage, bad_coverage, randomization_rate, nseconds, use_pbar):
    """
    @param input_lines: lines of input of csv data including the header
    @param good_coverage: the expected number of reads at informative positions
    @param bad_coverage: the expected number of reads at uninformative positions
    @param randomization_rate: the probability of an error per read
    @param nseconds: None or impose a time limit of this many seconds
    @param use_pbar: True iff a progress bar should be used
    @return: a multi-line string of the annotated csv file
    """
    verbose = False
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
    # define the three models
    homozygous = ReadCoverage.Homozygous(randomization_rate, good_coverage)
    heterozygous = ReadCoverage.Heterozygous(randomization_rate, good_coverage)
    overcovered = ReadCoverage.Overcovered(randomization_rate, bad_coverage)
    models = [homozygous, heterozygous, overcovered]
    # initialize the output header row
    header_row = input_rows[0]
    output_header_row = header_row[:5]
    for heading in header_row[5:]:
        if heading.endswith('sco'):
            output_header_row.append(heading)
        elif heading.endswith('cov'):
            output_header_row.extend([heading, heading+'_hom', heading+'_het', heading+'_ovr'])
        else:
            raise ValueError('each heading after the fifth should end with sco or cov')
    # get the rest of the rows
    data_rows = input_rows[1:]
    # define the number of genomic positions and the number of strains
    npositions = len(data_rows) / 5
    nstrains = (ncolumns - 5) / 2
    # begin the output
    out = StringIO()
    print >> out, ','.join(output_header_row)
    # initialize some stuff
    start_time = time.time()
    pbar = Progress.Bar(npositions) if use_pbar else None
    try:
        for position in range(npositions):
            # check the time
            if nseconds and time.time() - start_time > nseconds:
                raise TimeoutError()
            # get a chunk of five consecutive rows
            position_rows = [data_rows[position*5 + i] for i in range(5)]
            # get the corresponding log likelihoods
            log_likelihood_lists = get_log_likelihoods_per_strain(position_rows, models)
            # construct five annotated output lines
            for position_row in position_rows:
                output_row = position_row[:5]
                for i, log_likelihoods in enumerate(log_likelihood_lists):
                    # add the coverage, three annotations, and the score
                    coverage_string = position_row[5 + 2*i]
                    score_string = position_row[5 + 2*i + 1]
                    if log_likelihoods:
                        annotations = [str(x) for x in log_likelihoods]
                    else:
                        annotations = ['-', '-', '-']
                    output_row.extend([coverage_string] + annotations + [score_string])
                print >> out, ','.join(output_row)
            # update the progress bar
            if pbar:
                pbar.update(position + 1)
    except KeyboardInterrupt, e:
        if pbar:
            pbar.finish()
        raise e
    except TimeoutError, e:
        if pbar:
            pbar.finish()
        raise e
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
    print process(lines, options.good_coverage, options.bad_coverage, options.randomization_rate, nseconds, use_pbar)

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('--good-coverage', dest='good_coverage', type='int', default=20, help='expected read coverage of informative positions')
    parser.add_option('--bad-coverage', dest='bad_coverage', type='int', default=100, help='expected read coverage of overcovered positions')
    parser.add_option('--randomization-rate', dest='randomization_rate', type='float', default=0.1, help='randomization probability per read')
    options, args = parser.parse_args()
    main(options)

