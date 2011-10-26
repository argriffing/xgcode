"""Do a MAPP analysis using MAPP.jar.
"""

from StringIO import StringIO
import string
import tempfile
import subprocess
import os

import argparse

from SnippetUtil import HandlingError
import Form
import FormOut
import Util
import FelTree
import Newick
import NewickIO
import MatrixUtil
import HtmlTable
import MAPP
import Codon
import LeafWeights
import const

g_annotation_data = const.read('20100526a')
g_tree_data = const.read('20100526b')

# TODO de-hardcode this path to the jar file
g_mapp_jar_path = '/home/argriffi/mapp-analysis/MAPP.jar'

g_ordered_taxon_names = [
        'hg18',
        'rheMac2',
        'mm9',
        'canFam2',
        'loxAfr2',
        'monDom4',
        'ornAna1',
        'galGal3',
        'anoCar1',
        'xenTro2',
        'gasAcu1']

g_unordered_taxon_names = set(g_ordered_taxon_names)


def get_reverse_complement(s):
    trans = {
            'A' : 'T',
            'T' : 'A',
            'C' : 'G',
            'G' : 'C'}
    return ''.join(trans[nt] for nt in reversed(s))

def lines_to_annotated_snps(lines):
    """
    Yield Snp objects.
    @param lines: raw lines
    """
    for clump in gen_clumped_lines(lines):
        try:
            yield Snp(clump)
        except KnownBadSnpError as e:
            pass

def gen_clumped_lines(lines):
    """
    Yield clumps of contiguous lines with no intervening blank lines.
    Also comment lines are removed.
    @param lines: raw lines
    """
    clump = []
    for line in lines:
        line = line.strip()
        if line.startswith('#'):
            continue
        if line:
            clump.append(line)
        elif clump:
            yield clump
            clump = []
    if clump:
        yield clump


class SnpError(Exception): pass

class KnownBadSnpError(SnpError): pass

class Snp(object):

    def __init__(self, lines):
        # check the number of lines
        if len(lines) != 4:
            msg_a = 'expected 4 lines per annotated SNP '
            msg_b = 'but found %d' % len(lines)
            raise SnpError(msg_a + msg_b)
        # check for a known bad snp
        if lines[-1] == '?':
            raise KnownBadSnpError()
        # unpack the lines into member variables
        self.parse_first_line(lines[0])
        self.codon = lines[1].upper()
        self.within_codon_pos = int(lines[2])
        self.column = [x.upper() if x.isalpha() else None for x in lines[3]]
        # do some basic validation
        if not self.column[0]:
            msg = 'expected an aligned human amino acid for each SNP'
            raise SnpError(msg)
        if self.codon not in Codon.g_non_stop_codons:
            msg = 'expected a codon but found ' + self.codon
            raise SnpError(msg)
        if self.within_codon_pos not in (1, 2, 3):
            msg_a = 'expected the within-codon position '
            msg_b = 'to be either 1, 2, or 3, '
            msg_c = 'but found %d' % self.within_codon_pos
            raise SnpError(msg_a + msg_b + msg_c)
        # Assert that the major allele is actually in the codon
        # at the correct position, taking into account strand orientation.
        expected_nt = self.major_allele
        if self.orientation == '+':
            codon = self.codon
            observed_nt = codon[self.within_codon_pos - 1]
        else:
            codon = get_reverse_complement(self.codon)
            observed_nt = codon[2 - (self.within_codon_pos - 1)]
        if expected_nt != observed_nt:
            raise SnpError(
                    'in ' + self.variant_id + ' the major allele '
                    'is not in the codon at the correct position '
                    'even after taking into account the strand orientation.')
        # Assert that the human amino acid,
        # the first amino acid in the column,
        # is equal to the translation of the codon
        # taking into account the strand orientation.
        if self.orientation == '+':
            codon = self.codon
        else:
            codon = get_reverse_complement(self.codon)
        expected_aa = self.column[0]
        observed_aa = Codon.g_codon_to_aa_letter[codon]
        if expected_aa != observed_aa:
            raise SnpError(
                    'in ' + self.variant_id + ' the oriented codon '
                    'does not translate to '
                    'the human amino acid in the aligned column')
        # Get the mutant amino acid.
        if self.orientation == '+':
            c = list(self.codon)
            c[self.within_codon_pos - 1] = self.minor_allele
        else:
            c = list(get_reverse_complement(self.codon))
            c[2 - (self.within_codon_pos - 1)] = self.minor_allele
        codon = ''.join(c)
        self.mutant_aa = Codon.g_codon_to_aa_letter[codon]
        # Assert that the mutation is a missense mutation,
        # that is, that the wild and mutant amino acids are different.
        if self.mutant_aa == self.column[0]:
            raise SnpError(
                    'expected a missense mutation '
                    'but the wild type amino acid '
                    'is the same as the mutant amino acid')

    def parse_first_line(self, line):
        """
        @param line: a line of comma separated values
        """
        # break the line into unquoted elements
        ignore = string.whitespace + '"'
        v = [x.strip(ignore) for x in line.split(',')]
        # check the number of elements on the line
        if len(v) != 8:
            msg_a = 'expected 8 elements on the first line '
            msg_b = 'but found %d' % len(v)
            raise SnpError(msg_a + msg_b)
        # unpack the elements into member variables
        self.variant_id = v[0]
        self.chromosome_name = v[1]
        self.position = int(v[2])
        self.gene_id = v[3]
        self.gene_name = v[4]
        self.major_allele = v[5].upper()
        self.minor_allele = v[6].upper()
        self.orientation = v[7]
        # do some basic validation
        if self.major_allele not in 'ACGT':
            msg = 'major allele is invalid nucleotide: ' + self.major_allele
            raise SnpError(msg)
        if self.minor_allele not in 'ACGT':
            msg = 'minor allele is invalid nucleotide: ' + self.minor_allele
            raise SnpError(msg)
        if self.orientation not in '+-':
            msg_a = 'expected the orientation to be + or - '
            msg_b = 'but found ' + self.orientation
            raise SnpError(msg_a + msg_b)

    def get_simple_column(self):
        return ''.join(('-' if not x else x) for x in self.column)


def get_form():
    """
    @return: the body of a form
    """
    default_tree = NewickIO.parse(g_tree_data, FelTree.NewickTree)
    default_tree_string = NewickIO.get_narrow_newick_string(default_tree, 60)
    # define the list of form objects
    form_objects = [
            Form.MultiLine('tree', 'tree',
                default_tree_string),
            Form.MultiLine('annotation', 'SNP annotations',
                g_annotation_data)]
    return form_objects

def get_form_out():
    return FormOut.Report()

def aa_letter_to_aa_index(aa_letter):
    """
    @param aa_letter: an amino acid letter
    @return: None or an index between 0 and 19
    """
    for i, aa in enumerate(Codon.g_aa_letters):
        if aa == aa_letter:
            return i

def get_response_content(fs):
    # get the list of annotated snps
    snps = list(lines_to_annotated_snps(StringIO(fs.annotation)))
    # define the names of the temporary files
    temp_tree_filename = tempfile.mktemp(suffix='.tree')
    temp_fasta_filename = tempfile.mktemp(suffix='.fa')
    # write the temporary tree file
    with open(temp_tree_filename, 'w') as fout:
        print >> fout, fs.tree
    # write the temporary fasta file
    columns = [snp.get_simple_column() for snp in snps]
    sequences = [''.join(row) for row in zip(*columns)]
    with open(temp_fasta_filename, 'w') as fout:
        for sequence, taxon in zip(sequences, g_ordered_taxon_names):
            print >> fout, '>' + taxon
            print >> fout, sequence
    # call the mapp program
    cmd = [
            'gij',
            '-jar',
            g_mapp_jar_path,
            '-t',
            temp_tree_filename,
            '-f',
            temp_fasta_filename]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    mapp_out, mapp_err = p.communicate()
    # get the mapp output snp lines and skip the header
    lines = mapp_out.splitlines()[1:]
    lines = [line.rstrip('\n') for line in lines]
    lines = [line for line in lines if line]
    if len(snps) != len(lines):
        msg_a = 'found %d snps ' % len(snps)
        msg_b = 'but %d mapp data output lines' % len(lines)
        raise HandlingError(msg_a + msg_b)
    # write the output
    out = StringIO()
    for data_line, snp in zip(lines, snps):
        ignore = string.whitespace + "'" + '"'
        row = [x.strip(ignore) for x in data_line.split('\t')]
        if len(row) != 54:
            msg_a = 'expected 54 mapp output columns '
            msg_b = 'but found %d' % len(row)
            raise HandlingError(msg_a + msg_b + '\n' + str(row))
        impacts = [float(x) for x in row[12:32]]
        pvalues = [float(x) for x in row[32:52]]
        impact = impacts[aa_letter_to_aa_index(snp.mutant_aa)]
        pvalue = pvalues[aa_letter_to_aa_index(snp.mutant_aa)]
        trans = snp.column[0] + '->' + snp.mutant_aa
        col = snp.get_simple_column()
        row_out = [snp.variant_id, trans, col, impact, pvalue]
        print >> out, '\t'.join(str(x) for x in row_out)
    # return the response
    return out.getvalue()

def main(args):
    pass

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    args = parser.parse_args()
    main(args)
