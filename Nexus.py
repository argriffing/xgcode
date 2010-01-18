"""Read and write nexus files.
"""

import unittest
import StringIO
import Fasta
import Newick
import Monospace
import Util

nexus_sample_string = """
#NEXUS

[category 1: weight: 3, kappa: 2, A : 1, C : 1, G : 1, T : 1]
[category 2: weight: 1, kappa: 2, A : 1, C : 4, G : 4, T : 1]

BEGIN TAXA;
  DIMENSIONS ntax = 5;
  TAXLABELS Human Gorilla Gibbon Orangutan Chimpanzee;
END;

BEGIN TREES;
  TREE primates = (((Human:0.1, Chimpanzee:0.2):0.8, Gorilla:0.3):0.7, Orangutan:0.4, Gibbon:0.5);
END;

BEGIN CHARACTERS;
  DIMENSIONS nchar = 100;
  FORMAT datatype = DNA;
  MATRIX
    Human      CGGATGACGGTCTGGCGAATCCCCCTTCATGACCTTTGTTGGTCTAAGGTTGCGGGTTAGCCTAATGCACAGGCGGGGGCCGGTCGGGCGCCGCACCGTC
    Gorilla    CAAAACATTTGCGGAAACGTGAGCCCTCGTCTCCTGCGTCGCCTCAGTGAGCCGTCACCTCGACCAAGGCGGACGCCCGCCGGTTCGCCGCCGTGCCGCG
    Gibbon     GGCAGTATTTCCATTGCGATGGTTCACTTACACCGGTGAACACCAAATTCGCAGCCACGTCATCAAGCGGGCGTCTAAGCGGAGTGGTTGGCTGGCCCCC
    Orangutan  CGTTCGATCTGCTCCTCGGTAGTCCCCTTACACCTGTGGTACCCAAGTTTACGGCCCCCGCAAGCCGAGGTCGCCCACGTGGTTTCGACGTTCGGCGCTC
    Chimpanzee CGGATGACCGTCTGGCGGATCACTCTTCATGTGCTTTGTTAGTGTAAGGGTGCGCGTCAGCCATTCGCACGGGCAGGGACCGGTTGGGCGCCGCACCCAC
  ;
END;
"""


class NexusError(Exception):
    pass


class Nexus:
    """
    As far as I am concerned, this is a data structure that aggregates an alignment and a tree.
    """

    def __init__(self):
        # The nexus object should have an alignment and a tree.
        self.alignment = None
        self.tree = None
        # Keep an array of lines of comments.
        self.comments = []

    def add_comment(self, comment):
        """
        @param comment: a comment that is a string that is a single line
        """
        self.comments.append(comment)

    def load(self, lines):
        """
        @param lines: lines of nexus data
        """
        # get the taxa, tree, and character lines
        taxa_lines = []
        tree_lines = []
        character_lines = []
        current_array = None
        for line in Util.stripped_lines(lines):
            # Ignore an entire line that is a comment.
            # Nested comments and multi-line comments are not correctly processed here.
            if line.startswith('[') and line.endswith(']'):
                self.add_comment(line[1:-1])
                continue
            tokens = line.upper().split()
            if tokens == ['BEGIN', 'TAXA;']:
                current_array = taxa_lines
            elif tokens == ['BEGIN', 'TREES;']:
                current_array = tree_lines
            elif tokens == ['BEGIN', 'CHARACTERS;']:
                current_array = character_lines
            elif tokens == ['END;']:
                current_array = None
            elif current_array is not None:
                current_array.append(line)
        # assert that tree lines and character lines are present
        if not tree_lines:
            raise NexusError('TREES was not found')
        if not character_lines:
            raise NexusError('CHARACTERS was not found')
        # read the newick tree string
        nexus_tree_string = ''.join(tree_lines)
        if nexus_tree_string.count(';') != 1:
            raise NexusError('expected exactly one semicolon in the nexus TREES block')
        if nexus_tree_string.count('=') != 1:
            raise NexusError('expected exactly one equals sign in the nexus TREES block')
        offset = nexus_tree_string.find('=')
        newick_string = nexus_tree_string[offset+1:]
        self.tree = Newick.parse(newick_string, Newick.NewickTree)
        # read the alignment matrix
        arr = []
        found_matrix = False
        for line in character_lines:
            if line.upper().startswith('DIMENSIONS'):
                continue
            if line.upper().startswith('FORMAT'):
                continue
            if line.upper().startswith('MATRIX'):
                found_matrix = True
                continue
            if found_matrix:
                arr.append(line.replace(';', ' '))
        if not arr:
            raise NexusError('no alignment was found')
        tokens = ' '.join(arr).split()
        if len(tokens) % 2 != 0:
            raise NexusError('expected the alignment to be a list of (taxon, sequence) pairs')
        alignment_out = StringIO.StringIO()
        for header, sequence in Util.chopped(tokens, 2):
            sequence = sequence.upper()
            unexpected_letters = set(sequence) - set('ACGT')
            if unexpected_letters:
                raise NexusError('unexpected sequence character(s): %s' % list(unexpected_letters))
            print >> alignment_out, '>%s' % header
            print >> alignment_out, sequence
        alignment_string = alignment_out.getvalue()
        self.alignment = Fasta.Alignment(StringIO.StringIO(alignment_string))

    def __str__(self):
        out = StringIO.StringIO()
        alignment = self.alignment
        tree = self.tree
        # write the taxa block
        print >> out, '#NEXUS'
        print >> out, ''
        for comment in self.comments:
            print >> out, '[%s]' % comment
        print >> out, ''
        print >> out, 'BEGIN TAXA;'
        print >> out, '  DIMENSIONS ntax = %d;' % len(alignment.headers)
        print >> out, '  TAXLABELS %s;' % ' '.join(alignment.headers)
        print >> out, 'END;'
        # write the tree block
        print >> out, ''
        print >> out, 'BEGIN TREES;'
        print >> out, '  TREE primates = %s' % tree.get_newick_string()
        print >> out, 'END;'
        # write the alignment block
        print >> out, ''
        print >> out, 'BEGIN CHARACTERS;'
        print >> out, '  DIMENSIONS nchar = %d;' % len(alignment.columns)
        print >> out, '  FORMAT datatype = DNA;'
        print >> out, '  MATRIX'
        max_header_length = max(len(header) for header in alignment.headers)
        for header, sequence in zip(alignment.headers, alignment.sequences):
            print >> out, '    %s %s' % (Monospace.left_justify(header, max_header_length, ' '), sequence)
        print >> out, '  ;'
        print >> out, 'END;'
        return out.getvalue()


def get_sample_nexus_object():
    nexus = Nexus()
    nexus.load(StringIO.StringIO(nexus_sample_string))
    return nexus


class TestNexus(unittest.TestCase):

    def test_nexus(self):
        return get_sample_nexus_object()


def main():
    print get_sample_nexus_object()

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    #parser.add_option('-v', '--verbose', action='store_true', dest='verbose', default=False)
    #parser.add_option('-o', '--output', dest='output_filename', metavar='FILE', help='output file')
    parser.add_option('--test', action='store_true', dest='test', default=False)
    options, args = parser.parse_args()
    if options.test:
        suite = unittest.TestLoader().loadTestsFromTestCase(TestNexus)
        unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        main()
