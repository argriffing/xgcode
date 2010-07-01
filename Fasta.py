"""
Read and write FASTA files.
"""

from StringIO import StringIO
import unittest

import Codon
import Util
import iterutils


# from wikipedia
example_fasta_unaligned = """
>SEQUENCE_1
MTEITAAMVKELRESTGAGMMDCKNALSETNGDFDKAVQLLREKGLGKAAKKADRLAAEG
LVSVKVSDDFTIAAMRPSYLSYEDLDMTFVENEYKALVAELEKENEERRRLKDPNKPEHK
IPQFASRKQLSDAILKEAEEKIKEELKAQGKPEKIWDNIIPGKMNSFIADNSQLDSKLTL
MGQFYVMDDKKTVEQVIAEKEKEFGGKIKIVEFICFEVGEGLEKKTEDFAAEVAAQL
>SEQUENCE_2
SATVSEINSETDFVAKNDQFIALTKDTTAHIQSNSLQSVEELHSSTINGVKFEEYLKSQI
ATIGENLVVRRFATLKAGANGVVNGYIHTNGRVGVVIAAACDSAEVASKSRDLLRQICMH
"""

# from the bioperl wiki
example_fasta_aligned = """
>CYS1_DICDI
-----MKVILLFVLAVFTVFVSS---------------RGIPPEEQ------------SQ
FLEFQDKFNKKY-SHEEYLERFEIFKSNLGKIEELNLIAINHKADTKFGVNKFADLSSDE
FKNYYLNNKEAIFTDDLPVADYLDDEFINSIPTAFDWRTRG-AVTPVKNQGQCGSCWSFS
TTGNVEGQHFISQNKLVSLSEQNLVDCDHECMEYEGEEACDEGCNGGLQPNAYNYIIKNG
GIQTESSYPYTAETGTQCNFNSANIGAKISNFTMIP-KNETVMAGYIVSTGPLAIAADAV
E-WQFYIGGVF-DIPCN--PNSLDHGILIVGYSAKNTIFRKNMPYWIVKNSWGADWGEQG
YIYLRRGKNTCGVSNFVSTSII--
>ALEU_HORVU
MAHARVLLLALAVLATAAVAVASSSSFADSNPIRPVTDRAASTLESAVLGALGRTRHALR
FARFAVRYGKSYESAAEVRRRFRIFSESLEEVRSTN----RKGLPYRLGINRFSDMSWEE
FQATRL-GAAQTCSATLAGNHLMRDA--AALPETKDWREDG-IVSPVKNQAHCGSCWTFS
TTGALEAAYTQATGKNISLSEQQLVDCAGGFNNF--------GCNGGLPSQAFEYIKYNG
GIDTEESYPYKGVNGV-CHYKAENAAVQVLDSVNITLNAEDELKNAVGLVRPVSVAFQVI
DGFRQYKSGVYTSDHCGTTPDDVNHAVLAVGYGVENGV-----PYWLIKNSWGADWGDNG
YFKMEMGKNMCAIATCASYPVVAA
>CATH_HUMAN
------MWATLPLLCAGAWLLGV--------PVCGAAELSVNSLEK------------FH
FKSWMSKHRKTY-STEEYHHRLQTFASNWRKINAHN----NGNHTFKMALNQFSDMSFAE
IKHKYLWSEPQNCSAT--KSNYLRGT--GPYPPSVDWRKKGNFVSPVKNQGACGSCWTFS
TTGALESAIAIATGKMLSLAEQQLVDCAQDFNNY--------GCQGGLPSQAFEYILYNK
GIMGEDTYPYQGKDGY-CKFQPGKAIGFVKDVANITIYDEEAMVEAVALYNPVSFAFEVT
QDFMMYRTGIYSSTSCHKTPDKVNHAVLAVGYGEKNGI-----PYWIVKNSWGPQWGMNG
YFLIERGKNMCGLAACASYPIPLV
"""

# from brown.nuc in paml
brown_example_alignment = """
>Human
AAGCTTCACCGGCGCAGTCATTCTCATAATCGCCCACGGACTTACATCCTCATTACTATT
CTGCCTAGCAAACTCAAACTACGAACGCACTCACAGTCGCATCATAATCCTCTCTCAAGG
ACTTCAAACTCTACTCCCACTAATAGCTTTTTGATGACTTCTAGCAAGCCTCGCTAACCT
CGCCTTACCCCCCACTATTAACCTACTGGGAGAACTCTCTGTGCTAGTAACCACGTTCTC
CTGATCAAATATCACTCTCCTACTTACAGGACTCAACATACTAGTCACAGCCCTATACTC
CCTCTACATATTTACCACAACACAATGGGGCTCACTCACCCACCACATTAACAACATAAA
ACCCTCATTCACACGAGAAAACACCCTCATGTTCATACACCTATCCCCCATTCTCCTCCT
ATCCCTCAACCCCGACATCATTACCGGGTTTTCCTCTTGTAAATATAGTTTAACCAAAAC
ATCAGATTGTGAATCTGACAACAGAGGCTTACGACCCCTTATTTACCGAGAAAGCTCACA
AGAACTGCTAACTCATGCCCCATGTCTGACAACATGGCTTTCTCAACTTTTAAAGGATA
ACAGCTATCCATTGGTCTTAGGCCCCAAAAATTTTGGTGCAACTCCAAATAAAAGTAATA
ACCATGCACACTACTATAACCACCCTAACCCTGACTTCCCTAATTCCCCCCATCCTTACC
ACCCTCGTTAACCCTAACAAAAAAAACTCATACCCCCATTATGTAAAATCCATTGTCGCA
TCCACCTTTATTATCAGTCTCTTCCCCACAACAATATTCATGTGCCTAGACCAAGAAGTT
ATTATCTCGAACTGACACTGAGCCACAACCCAAACAACCCAGCTCTCCCTAAGCTT    
>Chimpanzee
AAGCTTCACCGGCGCAATTATCCTCATAATCGCCCACGGACTTACATCCTCATTATTATT
CTGCCTAGCAAACTCAAATTATGAACGCACCCACAGTCGCATCATAATTCTCTCCCAAGG
ACTTCAAACTCTACTCCCACTAATAGCCTTTTGATGACTCCTAGCAAGCCTCGCTAACCT
CGCCCTACCCCCTACCATTAATCTCCTAGGGGAACTCTCCGTGCTAGTAACCTCATTCTC
CTGATCAAATACCACTCTCCTACTCACAGGATTCAACATACTAATCACAGCCCTGTACTC
CCTCTACATGTTTACCACAACACAATGAGGCTCACTCACCCACCACATTAATAACATAAA
GCCCTCATTCACACGAGAAAATACTCTCATATTTTTACACCTATCCCCCATCCTCCTTCT
ATCCCTCAATCCTGATATCATCACTGGATTCACCTCCTGTAAATATAGTTTAACCAAAAC
ATCAGATTGTGAATCTGACAACAGAGGCTCACGACCCCTTATTTACCGAGAAAGCTTATA
AGAACTGCTAATTCATATCCCATGCCTAACAACATGGCTTTCTCAACTTTTAAAGGATA
ACAGCCATCCGTTGGTCTTAGGCCCCAAAAATTTTGGTGCAACTCCAAATAAAAGTAATA
ACCATGTATACTACCATAACCACCTTAACCCTAACTCCCTTAATTCTCCCCATCCTCACC
ACCCTCATTAACCCTAACAAAAAAAACTCATATCCCCATTATGTGAAATCCATTATCGCG
TCCACCTTTATCATTAGCCTTTTCCCCACAACAATATTCATATGCCTAGACCAAGAAGCT
ATTATCTCAAACTGGCACTGAGCAACAACCCAAACAACCCAGCTCTCCCTAAGCTT    
>Gorilla   
AAGCTTCACCGGCGCAGTTGTTCTTATAATTGCCCACGGACTTACATCATCATTATTATT
CTGCCTAGCAAACTCAAACTACGAACGAACCCACAGCCGCATCATAATTCTCTCTCAAGG
ACTCCAAACCCTACTCCCACTAATAGCCCTTTGATGACTTCTGGCAAGCCTCGCCAACCT
CGCCTTACCCCCCACCATTAACCTACTAGGAGAGCTCTCCGTACTAGTAACCACATTCTC
CTGATCAAACACCACCCTTTTACTTACAGGATCTAACATACTAATTACAGCCCTGTACTC
CCTTTATATATTTACCACAACACAATGAGGCCCACTCACACACCACATCACCAACATAAA
ACCCTCATTTACACGAGAAAACATCCTCATATTCATGCACCTATCCCCCATCCTCCTCCT
ATCCCTCAACCCCGATATTATCACCGGGTTCACCTCCTGTAAATATAGTTTAACCAAAAC
ATCAGATTGTGAATCTGATAACAGAGGCTCACAACCCCTTATTTACCGAGAAAGCTCGTA
AGAGCTGCTAACTCATACCCCGTGCTTAACAACATGGCTTTCTCAACTTTTAAAGGATA
ACAGCTATCCATTGGTCTTAGGACCCAAAAATTTTGGTGCAACTCCAAATAAAAGTAATA
ACTATGTACGCTACCATAACCACCTTAGCCCTAACTTCCTTAATTCCCCCTATCCTTACC
ACCTTCATCAATCCTAACAAAAAAAGCTCATACCCCCATTACGTAAAATCTATCGTCGCA
TCCACCTTTATCATCAGCCTCTTCCCCACAACAATATTTCTATGCCTAGACCAAGAAGCT
ATTATCTCAAGCTGACACTGAGCAACAACCCAAACAATTCAACTCTCCCTAAGCTT    
>Orangutan 
AAGCTTCACCGGCGCAACCACCCTCATGATTGCCCATGGACTCACATCCTCCCTACTGTT
CTGCCTAGCAAACTCAAACTACGAACGAACCCACAGCCGCATCATAATCCTCTCTCAAGG
CCTTCAAACTCTACTCCCCCTAATAGCCCTCTGATGACTTCTAGCAAGCCTCACTAACCT
TGCCCTACCACCCACCATCAACCTTCTAGGAGAACTCTCCGTACTAATAGCCATATTCTC
TTGATCTAACATCACCATCCTACTAACAGGACTCAACATACTAATCACAACCCTATACTC
TCTCTATATATTCACCACAACACAACGAGGTACACCCACACACCACATCAACAACATAAA
ACCTTCTTTCACACGCGAAAATACCCTCATGCTCATACACCTATCCCCCATCCTCCTCTT
ATCCCTCAACCCCAGCATCATCGCTGGGTTCGCCTACTGTAAATATAGTTTAACCAAAAC
ATTAGATTGTGAATCTAATAATAGGGCCCCACAACCCCTTATTTACCGAGAAAGCTCACA
AGAACTGCTAACTCTCACTCCATGTGTAACAACATGGCTTTCTCAGCTTTTAAAGGATA
ACAGCTATCCCTTGGTCTTAGGATCCAAAAATTTTGGTGCAACTCCAAATAAAAGTAACA
GCCATGTTTACCACCATAACTGCCCTCACCTTAACTTCCCTAATCCCCCCCATTACCGCT
ACCCTCATTAACCCCAACAAAAAAAACCCATACCCCCACTATGTAAAAACGGCCATCGCA
TCCGCCTTTACTATCAGCCTTATCCCAACAACAATATTTATCTGCCTAGGACAAGAAACC
ATCGTCACAAACTGATGCTGAACAACCACCCAGACACTACAACTCTCACTAAGCTT    
>Gibbon    
AAGCTTTACAGGTGCAACCGTCCTCATAATCGCCCACGGACTAACCTCTTCCCTGCTATT
CTGCCTTGCAAACTCAAACTACGAACGAACTCACAGCCGCATCATAATCCTATCTCGAGG
GCTCCAAGCCTTACTCCCACTGATAGCCTTCTGATGACTCGCAGCAAGCCTCGCTAACCT
CGCCCTACCCCCCACTATTAACCTCCTAGGTGAACTCTTCGTACTAATGGCCTCCTTCTC
CTGGGCAAACACTACTATTACACTCACCGGGCTCAACGTACTAATCACGGCCCTATACTC
CCTTTACATATTTATCATAACACAACGAGGCACACTTACACACCACATTAAAAACATAAA
ACCCTCACTCACACGAGAAAACATATTAATACTTATGCACCTCTTCCCCCTCCTCCTCCT
AACCCTCAACCCTAACATCATTACTGGCTTTACTCCCTGTAAACATAGTTTAATCAAAAC
ATTAGATTGTGAATCTAACAATAGAGGCTCGAAACCTCTTGCTTACCGAGAAAGCCCACA
AGAACTGCTAACTCACTACCCATGTATAACAACATGGCTTTCTCAACTTTTAAAGGATA
ACAGCTATCCATTGGTCTTAGGACCCAAAAATTTTGGTGCAACTCCAAATAAAAGTAATA
GCAATGTACACCACCATAGCCATTCTAACGCTAACCTCCCTAATTCCCCCCATTACAGCC
ACCCTTATTAACCCCAATAAAAAGAACTTATACCCGCACTACGTAAAAATGACCATTGCC
TCTACCTTTATAATCAGCCTATTTCCCACAATAATATTCATGTGCACAGACCAAGAAACC
ATTATTTCAAACTGACACTGAACTGCAACCCAAACGCTAGAACTCTCCCTAAGCTT    
"""

# from a simulation
simulated_codon_alignment = """
>Gibbon
TGTGAAGTGCTGCCTTTATTCTACTCCCCTTACCCTAGATTCATTCGTTTAGAGGGTGCC
CAGATCAGACGGTTCATAATAGGCCCGGAAGGCCACGGAGGACGATATCCACCCCTCTTC
CTATGCAATTCATGTCCTAGGGTCTCTTAT
>Orangutan
TGTGAAGTGGTGCCTCTATTCTACACCCCTTATCGTCGATTCATTAGTCTAGAAGGTGCC
CAAATTAGACGGTTCATAATAGGCCCGCGAGGGTGTGCAGGTCGATTTCCACCCCTATAT
TCGTGCAATTTCTCCCCTAGGGTCTCCTAT
>Gorilla
TATGATGTACCGCCTTTACTCTACATCCCACACGCGAAGCTCATTAGTCCAGAAGAGTCC
CAGAACAGAGTGGTCGTAAGGACCGCTGGAAGCCATGCAGGTAGATTCCCGCCCTTTACC
CCGTGCAATGTCCTCCCTATTGTTTCATAT
>Chimpanzee
TGCGAGAGGGCTCGTTTACTGTACATGCCTTGCTCTCAACCTACTAGTCCAGGGGGGCCC
TTTAACACATCTACCGTGATGAGCTCCGTATGGCATGCAGGTAGATTCCCGCTCTTAACC
TCTTGTAATGTCTTCCGTTTTGTTTCATAT
>Human
TGCGAGAGGGCTCCTTTACTGTACATGCCTTACTCTCAACCCATTAGTCCAGAGGGGCCC
CTTAGCAAACTTACCGTGATGAGCTCCGTTTGGCATGCAGGTAGATTCCCGCCCTTAACC
TCCTGTAATGTCTTCCGTTTCGTTTCATAT
"""


class AlignmentError(Exception):
    pass

class NucleotideAlignmentError(AlignmentError):
    pass

class CodonAlignmentError(AlignmentError):
    pass

def create_alignment(headers, sequences):
    """
    @param headers: a list of headers
    @param sequences: a list of sequences
    @return: an Alignment object
    """
    aln = []
    for header, sequence in zip(headers, sequences):
        aln.append('>' + header)
        aln.append(sequence)
    return Alignment(StringIO('\n'.join(aln)))


class Alignment:

    def __init__(self, lines):
        header_sequence_pairs = list(gen_header_sequence_pairs(lines))
        self.headers = [header for header, sequence in header_sequence_pairs]
        self.sequences = [sequence for header, sequence in header_sequence_pairs]
        for header in self.headers:
            if not header:
                raise AlignmentError('each sequence should have a header')
        if len(set(len(sequence) for sequence in self.sequences)) != 1:
            raise AlignmentError('not all sequences are the same length')
        self.columns = zip(*self.sequences)

    def get_column_count(self):
        return len(self.columns)

    def get_sequence_count(self):
        return len(self.sequences)

    def get_column_multiset(self):
        """
        @return: a dictionary mapping each different column to the number of times that column appears
        """
        column_multiset = {}
        for column in self.columns:
            col = tuple(column)
            column_multiset[col] = column_multiset.get(col, 0) + 1
        return column_multiset

    def get_unique_letters(self):
        return ''.join(list(sorted(set(''.join(self.sequences)))))

    def force_nucleotide(self):
        """
        Force the alignment to be a nucleotide alignment.
        Capitalize each 'acgt' letter.
        Remove columns with a non 'ACGT' letter.
        """
        next_columns = []
        for column in self.columns:
            next_column = [c.upper() for c in column if c in 'ACGTacgt']
            if len(next_column) == len(column):
                next_columns.append(next_column)
        if not next_columns:
            raise AlignmentError('no aligned nucleotide columns were found')
        self.columns = next_columns
        self.sequences = [''.join(sequence) for sequence in zip(*self.columns)]

    def get_fasta_sequence(self, header, columns=60):
        """
        @param columns: the maximum number of residues per line
        @return: a string representing the header and sequence in fasta format
        """
        header_to_sequence = dict(zip(self.headers, self.sequences))
        sequence = header_to_sequence[header]
        arr = []
        arr.append('>' + header)
        arr.append('\n'.join(iterutils.chopped(sequence, columns)))
        return '\n'.join(arr)

    def to_fasta_string(self, columns=60):
        """
        @param columns: the maximum number of residues per line
        @return: a string representing the whole multiple sequence alignment in fasta format
        """
        return '\n'.join(self.get_fasta_sequence(header, columns) for header in self.headers)


class CodonAlignment:

    def __init__(self, lines):
        """
        @param lines: lines of a fasta file
        """
        if not lines:
            raise AlignmentError('no fasta lines were provided')
        header_sequence_pairs = list(gen_header_sequence_pairs(lines))
        self.headers = [header for header, sequence in header_sequence_pairs]
        nucleotide_sequences = [seq.upper() for header, seq in header_sequence_pairs]
        if not nucleotide_sequences:
            raise AlignmentError('no nucleotide sequences were found')
        for header in self.headers:
            if not header:
                raise AlignmentError('each sequence should have a header')
        if len(set(len(sequence) for sequence in nucleotide_sequences)) != 1:
            raise AlignmentError('not all sequences are the same length')
        for seq in nucleotide_sequences:
            invalid_states = set(seq) - set('ACGT-')
            if invalid_states:
                example_invalid_state = invalid_states.pop()
                raise CodonAlignmentError('invalid nucleotide: %s' % example_invalid_state)
        nucleotide_columns = zip(*nucleotide_sequences)
        if len(nucleotide_columns) % 3 != 0:
            raise CodonAlignmentError('the number of aligned nucleotide columns should be a multiple of three')
        gappy_codon_sequences = [list(iterutils.chopped(seq, 3)) for seq in nucleotide_sequences]
        if not gappy_codon_sequences:
            raise CodonAlignmentError('no codon sequences were found')
        observed_gappy_codons = set(Util.flattened_nonrecursive(gappy_codon_sequences))
        valid_gappy_codons = set(list(Codon.g_non_stop_codons) + ['---'])
        invalid_gappy_codons = observed_gappy_codons - valid_gappy_codons
        if invalid_gappy_codons:
            example_invalid_codon = invalid_gappy_codons.pop()
            raise CodonAlignmentError('invalid codon: %s' % example_invalid_codon)
        self.columns = [col for col in zip(*gappy_codon_sequences) if '---' not in col]
        if not self.columns:
            raise CodonAlignmentError('no ungapped codon columns were found')
        self.sequences = zip(*self.columns)

    def get_column_multiset(self):
        """
        @return: a dictionary mapping each different column to the number of times that column appears
        """
        column_multiset = {}
        for column in self.columns:
            col = tuple(column)
            column_multiset[col] = column_multiset.get(col, 0) + 1
        return column_multiset

    def get_fasta_sequence(self, header, columns=60):
        """
        @param columns: the maximum number of residues per line
        @return: a string representing the header and sequence in fasta format
        """
        header_to_sequence = dict(zip(self.headers, self.sequences))
        sequence = header_to_sequence[header]
        arr = []
        arr.append('>' + header)
        arr.append('\n'.join(''.join(codons) for codons in iterutils.chopped(sequence, columns)))
        return '\n'.join(arr)

    def to_fasta_string(self, columns=60):
        """
        @param columns: the maximum number of residues per line
        @return: a string representing the whole multiple sequence alignment in fasta format
        """
        return '\n'.join(self.get_fasta_sequence(header, columns) for header in self.headers)


def gen_header_sequence_pairs(lines, max_sequences=None):
    """
    Yields (header, sequence) pairs.
    @param lines: a source of fasta lines
    @param max_sequences: A cap on the number of sequences yielded.
    """
    if max_sequences == 0:
        return
    current_header = None
    current_seq_arr = []
    yielded_sequence_count = 0
    for line in lines:
        line = line.strip()
        if line.startswith(';'):
            pass
        elif line.startswith('>'):
            if current_header is not None:
                yield current_header, ''.join(current_seq_arr)
                yielded_sequence_count += 1
                if max_sequences is not None:
                    if yielded_sequence_count >= max_sequences:
                        return
                current_seq_arr = []
            current_header = line[1:].strip()
        elif current_header is not None:
            current_seq_arr.append(line)
    if current_header is not None:
        yield current_header, ''.join(current_seq_arr)


class TestFasta(unittest.TestCase):

    def setUp(self):
        nucleotide_string = '>foo\nacgt-CATac--ACGT\n>bar\nacgtACGTacgtA-GT'
        self.simple_alignment = Alignment(StringIO(nucleotide_string))
        self.simple_alignment.force_nucleotide()

    def test_alignment(self):
        fin = StringIO(brown_example_alignment)
        alignment = Alignment(fin)
        self.assertEquals(len(alignment.headers), 5)
        self.assertEquals(len(alignment.sequences), 5)
        self.assertEquals(len(alignment.columns), 895)
        self.assertEquals(alignment.get_unique_letters(), 'ACGT')

    def test_nucleotide(self):
        aln = self.simple_alignment
        self.assertEquals(len(aln.headers), 2)
        self.assertEquals(len(aln.sequences), 2)
        self.assertEquals(len(aln.columns), 12)
        self.assertEquals(aln.get_unique_letters(), 'ACGT')

    def test_fasta_output(self):
        expected_alignment_string = '>foo\nACGTCATACAGT\n>bar\nACGTCGTACAGT'
        observed_alignment_string = self.simple_alignment.to_fasta_string()
        self.assertEquals(expected_alignment_string, observed_alignment_string)

    def test_codon(self):
        alignment = CodonAlignment(StringIO(simulated_codon_alignment))
        self.assertEquals(len(alignment.sequences), 5)
        self.assertEquals(len(alignment.columns), 50)


def main():
    fin = StringIO(brown_example_alignment)
    print 'example alignment header and sequence length pairs:'
    for header, sequence in gen_header_sequence_pairs(fin):
        print header
        print len(sequence)


if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    #parser.add_option('-v', '--verbose', action='store_true', dest='verbose', default=False)
    #parser.add_option('-o', '--output', dest='output_filename', metavar='FILE', help='output file')
    parser.add_option('--test', action='store_true', dest='test', default=False)
    options, args = parser.parse_args()
    if options.test:
        suite = unittest.TestLoader().loadTestsFromTestCase(TestFasta)
        unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        main()




