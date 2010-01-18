"""
KGEA is an abbreviation of knownGene.exonAA.

When header lines are mentioned in this module, they look something like the following.
>uc009xhd.1_hg18_1_1 342 0 0 chr1:247178160-247179185+
This is a line in a UCSC variant of the Fasta format.
"""

import unittest
import sys
import os
import os.path

import Progress


class KGEAError(Exception):
    pass


def get_piece_index(fin, position):
    """
    @param fin: an index file that is open for reading
    @param position: the query position within the chromosome
    @return: None or the fasta piece index that corresponds to the given position
    """
    for line in fin:
        stripped_line = line.strip()
        if stripped_line:
            first_index_string, last_index_string, piece_index_string = stripped_line.split()
            first_index = int(first_index_string)
            last_index = int(last_index_string)
            piece_index = int(piece_index_string)
            if first_index <= position <= last_index:
                return piece_index

def get_line_list(fin, chromosome_string, chromosome_position):
    """
    @param fin: a fasta file that is open for reading
    @param chromosome_string: the name of the chromosome
    @param chromosome_position: the query position within the chromosome
    @return: None or lines of an alignment that include the requested position
    """
    for line_list in gen_line_lists(fin):
        header_line = line_list[0]
        p = LocationParser(header_line)
        if p.chromosome != chromosome_string:
            continue
        if p.first_index <= chromosome_position <= p.last_index:
            return line_list

def gen_line_lists(fin):
    """
    Yield lists of stripped lines that are separated by blank lines in the input file.
    This is optimized for a very long file with short runs of non-blank lines.
    @param fin: an open file for reading
    """
    lines = []
    for line in fin:
        stripped_line = line.strip()
        if stripped_line:
            lines.append(stripped_line)
        else:
            if lines:
                yield lines
                lines = []
    if lines:
        yield lines

def header_line_to_taxon(header_line):
    """
    @param header_line: a line in a UCSC formatted fasta file
    @return: the name of the taxon
    """
    split_by_space = header_line.split()
    assert len(split_by_space) > 1
    element = split_by_space[0]
    split_by_period = element.split('.')
    assert len(split_by_period) == 2
    split_by_underscore = split_by_period[1].split('_')
    assert len(split_by_underscore) > 2
    taxon = split_by_underscore[1]
    return taxon

def piece_index_to_filename(piece_index, original_filename):
    """
    @param piece_index: the index of the piece of the original file
    @param original_filename: the name of the original fasta file
    @return: a filename that combines the index and the input filename specified in the constructor
    """
    assert piece_index < 10000
    piece_index_string = '%04d' % piece_index
    original_filename_segments = original_filename.split('.')
    assert len(original_filename_segments) > 1
    new_filename_segments = original_filename_segments[:-1] + [piece_index_string] + [original_filename_segments[-1]]
    return '.'.join(new_filename_segments)

def piece_filename_to_index(piece_filename):
    """
    @param piece_filename: a filename that contains the piece index
    @return: the piece index
    """
    # read the piece index from the piece pathname
    name_fragments = piece_filename.split('.')
    assert len(name_fragments) > 2
    piece_index_string = name_fragments[-2]
    piece_index = int(piece_index_string)
    return piece_index

def gen_elements(piece_pathname):
    """
    Each yielded element is a tuple of three values.
    The first value is a chromosome string.
    The second value is a first index.
    The third value is an last index.
    @param piece_pathname: the path to a fasta piece of an original huge fasta file
    """
    # read a bunch of alignments from the piece of the original huge fasta file
    fin = open(piece_pathname)
    for line_list in gen_line_lists(fin):
        # We only care about the first line in each list,
        # and each one looks something like this:
        # >uc009xhd.1_hg18_1_1 342 0 0 chr1:247178160-247179185+
        header_line = line_list[0]
        p = LocationParser(header_line)
        yield (p.chromosome, p.first_index, p.last_index)
    fin.close()


class LocationParser:

    def __init__(self, header_line):
        """
        Read the genomic location from the header line, and set some member variables accordingly.
        Read the chromosome as a string.
        Read the first and last indices as integers.
        Read the strand as a '+' or '-' single character string.
        @param header_line: a line in a UCSC formatted fasta file
        """
        # initialize the variables that will be filled
        self.chromosome = None
        self.first_index = None
        self.last_index = None
        self.strand = None
        self.length = None
        # assert that the taxon is hg18
        taxon = header_line_to_taxon(header_line)
        assert taxon == 'hg18', 'the location parser is broken for non-hg18 lines'
        # split the line by spaces, preparing to process the second and the last element
        split_by_space = header_line.split()
        assert len(split_by_space) > 2
        element = split_by_space[-1]
        length_string = split_by_space[1]
        # read the length
        self.length = int(length_string)
        # read the items of interest from the last element
        strand = element[-1]
        assert strand in '-+', strand
        location_string = element[:-1]
        split_by_colon = location_string.split(':')
        assert len(split_by_colon) == 2, header_line
        chromosome_string, coordinate_string = split_by_colon
        split_coordinates = coordinate_string.split('-')
        assert len(split_coordinates) == 2, header_line
        first_index_string, last_index_string = split_coordinates
        first_index = int(first_index_string)
        last_index = int(last_index_string)
        # set the variables
        self.chromosome = chromosome_string
        self.first_index = first_index
        self.last_index = last_index
        self.strand = strand



def gen_line_list_lists(line_list_generator, approx_flattened_nlines):
    """
    Yields lists of lists of stripped lines such that the total number of lines has approximately some length.
    @param line_list_generator: yields lists of lines
    @param approx_flattened_nlines: a target for the number of lines to yield
    """
    nlines = 0
    line_list_list = []
    for line_list in line_list_generator:
        line_list_list.append(line_list)
        nlines += len(line_list)
        if nlines >= approx_flattened_nlines:
            yield line_list_list
            nlines = 0
            line_list_list = []
    if line_list_list:
        yield line_list_list


class Splitter:
    
    def __init__(self, filename, pieces_directory):
        """
        @param filename: the filename of the huge fasta file to split
        @param pieces_directory: the name of the directory where the pieces will go
        """
        self.target_directory = pieces_directory
        self.filename = filename
        self.chromosome_strings = set()
        self.approx_lines_per_piece = 10000

    def run(self, verbose=False):
        """
        This reads the input file and writes pieces to a directory.
        @param verbose: True if we want to report our progress
        """
        # if we are reporting our progress then get the size of the file
        if verbose:
            nbytes_total = os.path.getsize(self.filename)
            pbar = Progress.Bar(nbytes_total)
            nbytes_current_approx = 0
        # process the file, possibly updating the progress bar
        fin = open(self.filename)
        piece_index = 0
        for line_list_lists in gen_line_list_lists(gen_line_lists(fin), self.approx_lines_per_piece):
            piece_filename = piece_index_to_filename(piece_index, self.filename)
            piece_path = os.path.join(self.target_directory, piece_filename)
            fout = open(piece_path, 'w')
            for line_list in line_list_lists:
                print >> fout, '\n'.join(line_list)
                print >> fout
            fout.close()
            piece_index += 1
            if verbose:
                nbytes = sum(sum(len(line) for line in line_list) for line_list in line_list_lists)
                nbytes_current_approx += nbytes
                pbar.update(nbytes_current_approx)
        fin.close()
        # possibly stop the progress bar
        if verbose:
            pbar.finish()

    def dry_run(self):
        """
        This just reads the input file.
        """
        nchunks = 0
        nlines = 0
        fin = open(self.filename)
        for lines in gen_line_lists(fin):
            self.process_header_line(lines[0])
            nlines += len(lines)
            nchunks += 1
            if not nchunks % 10000:
                print 'nlines:', nlines
                print 'nchunks:', nchunks
                print
        print 'total nlines:', nlines
        print 'total nchunks:', nchunks
        fin.close()

    def process_header_line(self, header_line):
        coordinate_string = header_line.split()[-1].strip()
        chromosome_string = coordinate_string.split(':')[0]
        if chromosome_string not in self.chromosome_strings:
            self.chromosome_strings.add(chromosome_string)
            print 'new chromosome string in this header line:'
            print header_line
            print 'all chromosome strings:'
            for s in sorted(self.chromosome_strings):
                print s
            print 


class Indexer:
    """
    Index the fasta files.

    This script starts from the assumption that there is a subdirectory
    with a bunch of pieces of an original huge fasta file.
    The original fasta file was
    knownGene.exonAA.fa
    The pieces are in the directory
    fasta-pieces
    and are given names like
    knownGene.exonAA.1234.fa

    The job of this script is to index these fasta files by genomic coordinates.
    The hard parts are already done;
    in particular, the fasta header for each exon gives a chromosome name
    and the coordinates covered by the exon within the chromosome.

    This script simply makes, for each chromosome,
    a file that has three values in each row: the first coordinate, the last coordinate, and the piece index.
    Oh, it also makes a file that lists the valid chromosome names.
    """

    def __init__(self, index_directory, chromosome_list_filename, pieces_directory):
        """
        Initialize some parameters.
        @param index_directory: write index files to this directory
        @param chromosome_list_filename: write the chromosome strings to a file with this name
        @param pieces_directory: read the fasta files that are in this directory
        """
        self.index_directory = index_directory
        self.chromosome_list_filename = chromosome_list_filename
        self.pieces_directory = pieces_directory
        # check some conditions so we fail early instead of in the middle of a long run
        if not os.path.isdir(self.pieces_directory):
            raise KGEAError('missing the directory with fasta files: ' + self.pieces_directory)
        if not os.path.isdir(self.index_directory):
            raise KGEAError('missing the directory to which the index files should be written: ' + self.index_directory)

    def run(self, verbose=False):
        """
        Create the index files.
        This might take a while.
        @param verbose: True if we want to write our progress to stdout
        """
        # fill a dictionary by reading all of the fasta pieces
        chromosome_string_to_rows = {}
        piece_filenames = list(sorted(os.listdir(self.pieces_directory)))
        if verbose:
            print >> sys.stderr, 'creating a dictionary from the fasta pieces:'
            pbar = Progress.Bar(len(piece_filenames))
            nfiles_read = 0
        for piece_filename in piece_filenames:
            piece_index = piece_filename_to_index(piece_filename)
            piece_pathname = os.path.join(self.pieces_directory, piece_filename)
            for chromosome_string, first_index, last_index in gen_elements(piece_pathname):
                row = (first_index, last_index, piece_index)
                rows = chromosome_string_to_rows.get(chromosome_string, [])
                rows.append(row)
                chromosome_string_to_rows[chromosome_string] = rows
            if verbose:
                nfiles_read += 1
                pbar.update(nfiles_read)
        # define the list of chromosome strings
        chromosome_strings = list(sorted(chromosome_string_to_rows))
        assert len(chromosome_strings) < 1000
        # write the list of valid chromosome strings
        fout = open(self.chromosome_list_filename, 'w')
        fout.write('\n'.join(chromosome_strings))
        fout.close()
        if verbose:
            print >> sys.stderr, 'wrote', self.chromosome_list_filename
            print >> sys.stderr, 'writing the index files:'
            pbar = Progress.Bar(len(chromosome_strings))
            nwritten = 0
        # for each chromosome string write the index file
        for chromosome_string in chromosome_strings:
            rows = chromosome_string_to_rows[chromosome_string]
            index_filename = chromosome_string + '.index'
            index_pathname = os.path.join(self.index_directory, index_filename)
            fout = open(index_pathname, 'w')
            for row in sorted(rows):
                print >> fout, '%d\t%d\t%d' % row
            fout.close()
            if verbose:
                nwritten += 1
                pbar.update(nwritten)


class Finder:

    def __init__(self, index_directory, chromosome_list_filename, fasta_directory):
        """
        Initialize some parameters.
        @param index_directory: index files have been written to this directory
        @param chromosome_list_filename: the chromosome strings have been written to a file with this name
        @param fasta_directory: the fasta files are in this directory
        """
        self.index_directory = index_directory
        self.fasta_directory = fasta_directory
        self.valid_chromosome_strings = None
        # check some conditions so we fail early
        if not os.path.isdir(self.fasta_directory):
            raise KGEAError('missing the directory with fasta files: ' + self.fasta_directory)
        if not os.path.isdir(self.index_directory):
            raise KGEAError('missing the directory to which the index files should be written: ' + self.index_directory)
        # try to get the list of valid chromosome strings
        try:
            fin = open(chromosome_list_filename)
        except IOError, e:
            raise KGEAError('there was a problem opening the chromosome list file: ' + chromosome_list_filename)
        lines = [line.strip() for line in fin.readlines()]
        fin.close()
        self.valid_chromosome_strings = set(line for line in lines if line)

    def get_alignment_lines(self, chromosome_string, chromosome_position, verbose=False):
        """
        @param chromosome_string: a string like 'chr7' or 'chrX'
        @param chromosome_position: a non-negative integer nucleotide offset
        @param verbose: True if we want to write our progress to stdout
        @return: a list of lines of the fasta alignment of the coding part of some exon
        """
        # validate the chromosome string
        if chromosome_string not in self.valid_chromosome_strings:
            raise KGEAError('invalid chromosome: ' + chromosome_string)
        # validate the chromosome position
        if chromosome_position < 0:
            raise KGEAError('the nucleotide offset on the chromosome cannot be negative')
        # get the index pathname using the chromosome string
        index_filename = chromosome_string + '.index'
        index_pathname = os.path.join(self.index_directory, index_filename)
        # read the index file to find a piece index
        if verbose:
            print 'searching the index file', index_pathname
        fin = open(index_pathname)
        piece_index = get_piece_index(fin, chromosome_position)
        fin.close()
        if not piece_index:
            return []
        # read the fasta piece corresponding to the index
        piece_pathname = self.fasta_directory + '/knownGene.exonAA.%04d.fa' % piece_index
        if verbose:
            print 'searching the fasta file', piece_pathname
        fin = open(piece_pathname)
        line_list = get_line_list(fin, chromosome_string, chromosome_position)
        fin.close()
        if line_list is None:
            raise KGEAError('the index pointed to a fasta file that does not have the requested position')
        # return the lines of the alignment
        return line_list

    def gen_taxon_aa_pairs(self, chromosome_string, chromosome_position, verbose=False):
        """
        Yield (taxon name, amino acid) pairs
        @param chromosome_string: a string like 'chr7' or 'chrX'
        @param chromosome_position: a non-negative integer nucleotide offset
        @param verbose: True if we want to write our progress to stdout
        """
        # get the alignment list associated with the position on the chromosome
        fasta_line_list = self.get_alignment_lines(chromosome_string, chromosome_position, verbose)
        if not fasta_line_list:
            return
        ntaxa, remainder = divmod(len(fasta_line_list), 2)
        if remainder:
            raise KGEAError('the number of lines in each alignment was expected to be even')
        # get the amino acid offset of interest
        header_line = fasta_line_list[0]
        p = LocationParser(header_line)
        aa_offset = (chromosome_position - p.first_index) / 3
        # yield lines at the column of interest
        nspecies = len(fasta_line_list) / 2
        for i in range(nspecies):
            # get the species name
            header_line = fasta_line_list[i*2]
            taxon_name = header_line_to_taxon(header_line)
            # get the amino acid
            if p.strand == '+':
                sequence_line = fasta_line_list[i*2+1]
            elif p.strand == '-':
                sequence_line = ''.join(reversed(fasta_line_list[i*2+1]))
            aa = sequence_line[aa_offset]
            # show this row of the alignment column
            yield (taxon_name, aa)

    def get_column_lines(self, chromosome_string, chromosome_position, verbose=False):
        """
        @param chromosome_string: a string like 'chr7' or 'chrX'
        @param chromosome_position: a non-negative integer nucleotide offset
        @param verbose: True if we want to write our progress to stdout
        @return: a list of column lines
        """
        return [taxon + '\t' + aa for taxon, aa in self.gen_taxon_aa_pairs(chromosome_string, chromosome_position, verbose)]


class TestParser(unittest.TestCase):

    def test_human_header(self):
        header_line = '>uc009xhd.1_hg18_1_1 342 0 0 chr1:247178160-247179185+'
        data = LocationParser(header_line)
        self.assertEqual(data.strand, '+')
        self.assertEqual(data.first_index, 247178160)
        self.assertEqual(data.last_index, 247179185)
        self.assertEqual(data.chromosome, 'chr1')
        taxon = header_line_to_taxon(header_line)
        self.assertEqual(taxon, 'hg18')

    def test_nonhuman_header(self):
        header_line = '>uc002jlv.1_gasAcu1_1_1 172 0 0 chrXI:16540840-16541230-'
        taxon = header_line_to_taxon(header_line)
        self.assertEqual(taxon, 'gasAcu1')


class TestSplitter(unittest.TestCase):
    pass


class TestIndexer(unittest.TestCase):
    pass


if __name__ == '__main__':
    suite = unittest.TestSuite([
        unittest.TestLoader().loadTestsFromTestCase(TestParser),
        unittest.TestLoader().loadTestsFromTestCase(TestSplitter),
        unittest.TestLoader().loadTestsFromTestCase(TestIndexer)])
    unittest.TextTestRunner(verbosity=2).run(suite)

