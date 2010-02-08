"""Break a resequencing file into files for individual strains.

This is not so useful as a internet accessible script.
Columns of each output file are defined as follows.
strain: a string defining the genetic line
chromosome: a string defining the name of the chromosome
offset: an integer offset in units of base pairs
A: the number of A reads aligned to the this position
C: the number of C reads aligned to the this position
G: the number of G reads aligned to the this position
T: the number of T reads aligned to the this position
gap: the number of gaps aligned to the this position
"""

from StringIO import StringIO
import time
import optparse
import sys
import os

import numpy as np

from SnippetUtil import HandlingError
import Form
import Progress
import iterutils


class TimeoutError(Exception): pass


g_output_header_row = [
        'genetic_line', 'chromosome', 'position',
        'A_count', 'C_count', 'G_count', 'T_count', 'gap_count']


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
    """

    def __init__(self, strain_name, chromosome_name, output_directory):
        """
        Start a new file.
        @param strain_name: the name of the associated genetic strain
        @param chromosome_name: the name of the chromosome
        @param output_directory: the directory where the parsed files will go
        """
        # save the strain and chromosome names
        self.strain_name = strain_name
        self.chromosome_name = chromosome_name
        # create a new open file for writing
        filename = chromosome_name + '_' + strain_name + '.csv'
        pathname = os.path.join(output_directory, filename)
        self.fout = open(pathname, 'w')
        # write the header
        print >> self.fout, ','.join(g_output_header_row)

    def finish(self):
        self.fout.close()

    def add_position(self, offset, coverages):
        """
        Add a position to the chromosome.
        @param offset: the nucleotide offset
        @param coverages: the (A, C, G, T, gap) coverage counts
        """
        # do validation
        if len(coverages) != 5:
            raise ValueError('five coverage counts were expected')
        # write a line to the open file
        row = [self.strain_name, self.chromosome_name, str(offset)]
        row.extend([str(x) for x in coverages])
        line = ','.join(row)
        print >> self.fout, line

def gen_validated_data_rows(row_source, header_row):
    """
    @param row_source: a source of non-header rows from a csv file
    @param header_row: the header row from a csv file
    """
    for row in row_source:
        if len(row) != len(header_row):
            raise ValueError('each row should have the same number of elements as the first row')
        yield row

def parse(row_source, output_directory_path):
    """
    @param row_source: a source of rows from a csv file
    @param output_directory_path: a place to put the output files
    @return: a list of chromosome objects
    """
    header_row = next(row_source)
    # validate the header row
    ncolumns = len(header_row)
    if ncolumns < 7:
        raise ValueError('there should be at least seven columns of input')
    if ncolumns % 2 != 1:
        raise ValueError('the number of input columns should be odd')
    # read the strain names from the header
    strain_names = []
    for heading in header_row[5:]:
        if heading.endswith('sco'):
            pass
        elif heading.endswith('cov'):
            strain_names.append(heading[:-3])
        else:
            raise ValueError('each heading after the fifth should end with sco or cov')
    # define a row generator that validates each row that it yields
    validated_row_source = gen_validated_data_rows(row_source, header_row)
    # go through validated rows in groups of five
    chromosome_dict = {}
    for position_rows in iterutils.grouper(validated_row_source, 5):
        if not all(position_rows):
            raise ValueError('the number of non-header rows was not a multiple of five')
        process_genomic_position(position_rows, strain_names, chromosome_dict)
    # turn the dictionary of chromosomes into a somewhat ordered list
    chromosomes = [chromosome for identifier, chromosome in sorted(chromosome_dict.items())]
    # return the list of chromosomes
    return chromosomes

def process_genomic_position(rows, strain_names, chromosome_dict):
    """
    The input data represents a single genomic position for multiple strains.
    Each chromosome object is accessible by its (strain_name, chromosome_name) pair.
    @param rows: five rows of the csv file
    @param strain_names: an ordered list of strain names
    @param chromosome_dict: a dictionary of chromosome objects
    """
    nheadercols = 5
    # Read the chromosome name from the first of the five rows.
    chromosome_name = rows[0][0]
    # Read the genomic position from the first of the five rows.
    genomic_position = int(rows[0][3])
    # Extract the important columns of the important rows,
    # and convert the elements from strings to integers.
    reduced_rows = []
    for row in rows:
        reduced_row = []
        for i, element in enumerate(row):
            if i >= nheadercols and i % 2 == 1:
                count = int(element)
                if count < 0:
                    raise ValueError('counts must be nonnegative: ' + str(count))
                reduced_row.append(count)
        reduced_rows.append(reduced_row)
    # Each reduced row corresponds to a nucleotide.
    # Each column in a reduced row corresponds to a strain.
    # Each element is a coverage count.
    coverages = zip(*reduced_rows)
    # Add an observation to each chromosome.
    for coverage, strain_name in zip(coverages, strain_names):
        # The input observation order is ACTG,
        # but the order I want to use is ACGT.
        A, C, T, G, gap = coverage
        coverage = (A, C, G, T, gap)
        # get or create the chromosome object
        identifier = (strain_name, chromosome_name)
        if identifier in chromosome_dict:
            chromosome = chromosome_dict[identifier]
        else:
            chromosome = Chromosome(strain_name, chromosome_name, output_directory_path)
            chromosome_dict[identifier] = chromosome
        # add the position to the chromosome
        chromosome.add_position(genomic_position, coverage)

def gen_rows(lines):
    """
    Yield rows from a source of csv lines.
    @param lines: an iterable source of lines
    """
    for line in lines:
        line = line.strip()
        if line:
            row = [x.strip() for x in line.split(',')]
            yield row

if __name__ == '__main__':
    from optparse import OptionParser
    usage = 'Usage: %prog <source.csv> <output-directory>'
    parser = OptionParser(usage=usage)
    options, args = parser.parse_args()
    input_filename, output_directory_path = args
    with open(input_filename) as fin:
        for chromosome in parse(gen_rows(fin), output_directory_path):
            chromosome.finish()
