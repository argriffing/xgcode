"""Create phylip files.
"""

import unittest
import StringIO

import Util
import Monospace
import Fasta


def get_alignment_string_non_interleaved(alignment):
    """
    @param alignment: a fasta alignment object
    @return: a non interleaved phylip alignment string
    """
    out = StringIO.StringIO()
    # print the number of sequences and the length of each sequence in the alignment
    print >> out, len(alignment.headers), len(alignment.sequences[0])
    # print each sequence
    for header, sequence in zip(alignment.headers, alignment.sequences):
        print >> out, header
        for segment in Util.chopped(sequence, 60):
            print >> out, segment
    return out.getvalue().strip()

def get_alignment_string_interleaved(alignment):
    """
    @param alignment: a fasta alignment object
    @return: an interleaved phylip alignment string
    """
    chopped_sequences = [Util.chopped(seq, 60) for seq in alignment.sequences]
    bands = zip(*chopped_sequences)
    out = StringIO.StringIO()
    lengths = [9] + [len(header) for header in alignment.headers]
    print >> out, len(alignment.headers), len(alignment.sequences[0])
    for header, segment in zip(alignment.headers, bands[0]):
        aligned_header = Monospace.left_justify(header, max(lengths), ' ')
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
    return Fasta.Alignment(StringIO.StringIO(Fasta.example_fasta_aligned))


class TestPhylip(unittest.TestCase):

    def test_foo(self):
        alignment = make_sample_alignment()
        get_alignment_string_interleaved(alignment)
        get_alignment_string_non_interleaved(alignment)


def main():
    alignment = make_sample_alignment()
    print get_alignment_string_non_interleaved(alignment)

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    #parser.add_option('-v', '--verbose', action='store_true', dest='verbose', default=False)
    #parser.add_option('-o', '--output', dest='output_filename', metavar='FILE', default='-', help='output file')
    parser.add_option('--test', action='store_true', dest='test', default=False, help='run some unit tests')
    options, args = parser.parse_args()
    # run a test or run a demo
    if options.test:
        suite = unittest.TestLoader().loadTestsFromTestCase(TestPhylip)
        unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        main()

