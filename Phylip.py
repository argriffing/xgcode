"""Create phylip files.

http://www.molecularevolution.org/resources/fileformats/phylip_dna
"""

import unittest
from StringIO import StringIO

import Fasta
import iterutils
import Util


class PhylipError(Exception): pass

def read_non_interleaved(raw_lines):
    """
    @param raw_lines: raw lines of a non-interleaved phylip alignment file
    @return: headers, sequences
    """
    lines = Util.get_stripped_lines(raw_lines)
    header_line, data_lines = lines[0], lines[1:]
    header_row = header_line.split()
    if len(header_row) != 2:
        raise PhylipError('the header should be a line with two integers')
    ntaxa_s, ncolumns_s = header_row
    try:
        ntaxa = int(ntaxa_s)
        ncolumns = int(ncolumns_s)
    except ValueError:
        raise PhylipError('the header should be a line with two integers')
    ntaxa_observed = len(data_lines)
    if ntaxa_observed != ntaxa:
        msg_a = 'the header says there are %d taxa' % ntaxa
        msg_b = 'but %d taxa were observed' % ntaxa_observed
        raise PhylipError(msg_a + msg_b)
    compound_data_rows = [x.split() for x in data_lines]
    for row in compound_data_rows:
        if len(row) != 2:
            raise PhylipError(
                    'expected a taxon and a sequence for each data row')
    headers, sequences = zip(*compound_data_rows)
    ncolumns_observed = len(sequences[0])
    for row in sequences:
        if len(row) != ncolumns_observed:
            raise PhylipError(
                    'expected the same number of columns on each row')
    if ncolumns_observed != ncolumns:
        msg_a = 'the header says there are %d alignment columns' % ncolumns
        msg_b = 'but %d alignment columns were observed' % ncolumns_observed
        raise PhylipError(msg_a + msg_b)
    return headers, sequences


def get_alignment_string_non_interleaved(alignment):
    """
    @param alignment: a fasta alignment object
    @return: a non interleaved phylip alignment string
    """
    out = StringIO()
    # print the number of sequences and the length of each sequence
    print >> out, len(alignment.headers), len(alignment.sequences[0])
    # print each sequence
    for header, sequence in zip(alignment.headers, alignment.sequences):
        print >> out, header
        for segment in iterutils.chopped(sequence, 60):
            print >> out, segment
    return out.getvalue().strip()

def get_alignment_string_interleaved(alignment):
    """
    @param alignment: a fasta alignment object
    @return: an interleaved phylip alignment string
    """
    chopped_sequences = [
            iterutils.chopped(seq, 60) for seq in alignment.sequences]
    bands = zip(*chopped_sequences)
    out = StringIO()
    print >> out, len(alignment.headers), len(alignment.sequences[0])
    lengths = [9] + [len(header) for header in alignment.headers]
    n = max(lengths)
    for header, segment in zip(alignment.headers, bands[0]):
        aligned_header = header.ljust(n)
        print >> out, '%s  %s' % (aligned_header, segment)
    print >> out, ''
    for band in bands[1:]:
        for segment in band:
            print >> out, segment
        print >> out, ''
    return out.getvalue().strip()

def make_sample_alignment():
    """
    Make a sample alignment object.
    """
    return Fasta.Alignment(StringIO(Fasta.example_fasta_aligned))


class TestPhylip(unittest.TestCase):

    def test_foo(self):
        alignment = make_sample_alignment()
        get_alignment_string_interleaved(alignment)
        get_alignment_string_non_interleaved(alignment)


if __name__ == '__main__':
    unittest.main()
