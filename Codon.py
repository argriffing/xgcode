"""This module provides names and relationships among amino acids, nucleotides, and codons.

For consistency all residue symbols, abbreviations, and names use capital letters.
"""

import unittest

class CodonError(Exception):
    pass

# begin external variable list

g_nt_letters = ()
g_nt_names = ()
g_aa_letters = ()
g_aa_abbrs = ()
g_aa_names = ()
g_all_codons = ()
g_stop_codons = ()
g_start_codon = 'ATG'
g_non_stop_codons = ()

g_sorted_non_stop_codons = ()
g_sorted_aa_letters = ()
g_sorted_nt_letters = ()

g_aa_name_to_aa_letter = {}
g_aa_abbr_to_aa_letter = {}
g_aa_letter_to_aa_name = {}
g_aa_letter_to_aa_abbr = {}
g_nt_letter_to_nt_name = {}
g_nt_name_to_nt_letter = {}

g_codon_to_aa_letter = {}
g_aa_letter_to_codons = {}
g_aa_letter_pair_to_min_codon_difference = {}
g_codon_to_missense_codons = {}

g_aa_pair_to_snp_count = {}

# end external variable list


_aa_codons = {
        'Ala': ('GCU', 'GCC', 'GCA', 'GCG'),
        'Leu': ('UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'),
        'Arg': ('CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'),
        'Lys': ('AAA', 'AAG'),
        'Asn': ('AAU', 'AAC'),     
        'Met': ('AUG',),
        'Asp': ('GAU', 'GAC'),  
        'Phe': ('UUU', 'UUC'),
        'Cys': ('UGU', 'UGC'), 
        'Pro': ('CCU', 'CCC', 'CCA', 'CCG'),
        'Gln': ('CAA', 'CAG'),
        'Ser': ('UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'),
        'Glu': ('GAA', 'GAG'),
        'Thr': ('ACU', 'ACC', 'ACA', 'ACG'),
        'Gly': ('GGU', 'GGC', 'GGA', 'GGG'),
        'Trp': ('UGG',),
        'His': ('CAU', 'CAC'),
        'Tyr': ('UAU', 'UAC'),
        'Ile': ('AUU', 'AUC', 'AUA'), 
        'Val': ('GUU', 'GUC', 'GUA', 'GUG')
        }

_aa_aliases = (
        ('A', 'ALA', 'alanine'),
        ('R', 'ARG', 'arginine'),
        ('N', 'ASN', 'asparagine'),
        ('D', 'ASP', 'aspartic acid'),
        ('C', 'CYS', 'cysteine'),
        ('Q', 'GLN', 'glutamine'),
        ('E', 'GLU', 'glutamic acid'),
        ('G', 'GLY', 'glycine'),
        ('H', 'HIS', 'histidine'),
        ('I', 'ILE', 'isoleucine'),
        ('L', 'LEU', 'leucine'),
        ('K', 'LYS', 'lysine'),
        ('M', 'MET', 'methionine'),
        ('F', 'PHE', 'phenylalanine'),
        ('P', 'PRO', 'proline'),
        ('S', 'SER', 'serine'),
        ('T', 'THR', 'threonine'),
        ('W', 'TRP', 'tryptophan'),
        ('Y', 'TYR', 'tyrosine'),
        ('V', 'VAL', 'valine')
        )

g_nt_letters = tuple('ACGT')
g_nt_names = ('ADENINE', 'CYTOSINE', 'GUANINE', 'THYMINE')
g_aa_letters = tuple(sorted(aa for aa, abbr, name in _aa_aliases))
g_aa_abbrs = tuple(abbr.upper() for aa, abbr, name in sorted(_aa_aliases))
g_aa_names = tuple(name.upper() for aa, abbr, name in sorted(_aa_aliases))
g_all_codons = tuple(a+b+c for a in g_nt_letters for b in g_nt_letters for c in g_nt_letters)
g_stop_codons = ('TAA', 'TGA', 'TAG')
g_non_stop_codons = tuple(codon for codon in g_all_codons if codon not in g_stop_codons)

g_sorted_non_stop_codons = tuple(sorted(g_non_stop_codons))
g_sorted_aa_letters = tuple(sorted(g_aa_letters))
g_sorted_nt_letters = tuple(sorted(g_nt_letters))

g_nt_letter_to_nt_name = dict((nt, name) for nt, name in zip(g_nt_letters, g_nt_names))
g_nt_name_to_nt_letter = dict((name, nt) for nt, name in zip(g_nt_letters, g_nt_names))
g_aa_letter_to_aa_name = dict((aa, name) for aa, name in zip(g_aa_letters, g_aa_names))
g_aa_letter_to_aa_abbr = dict((aa, abbr) for aa, abbr in zip(g_aa_letters, g_aa_abbrs))
g_aa_abbr_to_aa_letter = dict((abbr, aa) for aa, abbr in zip(g_aa_letters, g_aa_abbrs))
g_aa_name_to_aa_letter = dict((name, aa) for aa, name in zip(g_aa_letters, g_aa_names))

def _gen_aa_letter_and_codons():
    for aa_abbr, rna_codons in _aa_codons.items():
        aa_letter = g_aa_abbr_to_aa_letter[aa_abbr.upper()]
        cdna_codons = tuple(codon.replace('U', 'T') for codon in rna_codons)
        yield aa_letter, cdna_codons

def _gen_codon_and_aa_letter():
    for aa_letter, cdna_codons in _gen_aa_letter_and_codons():
        for codon in cdna_codons:
            yield codon, aa_letter

def _hamming_distance(first, second):
    return sum(1 for a, b in zip(first, second) if a != b)

def _gen_aa_letter_pair_and_min_codon_difference():
    for aa1 in g_aa_letters:
        c1 = g_aa_letter_to_codons[aa1]
        for aa2 in g_aa_letters:
            c2 = g_aa_letter_to_codons[aa2]
            yield (aa1, aa2), min(_hamming_distance(a, b) for a in c1 for b in c2)

def _gen_codon_and_missense_codons():
    for c1 in g_non_stop_codons:
        aa_letter = g_codon_to_aa_letter[c1]
        missense_list = []
        for c2 in g_non_stop_codons:
            if _hamming_distance(c1, c2) == 1:
                if g_codon_to_aa_letter[c2] != aa_letter:
                    missense_list.append(c2)
        yield c1, tuple(missense_list)

g_aa_letter_to_codons = dict(_gen_aa_letter_and_codons())
g_codon_to_aa_letter = dict(_gen_codon_and_aa_letter())
g_aa_letter_pair_to_min_codon_difference = dict(_gen_aa_letter_pair_and_min_codon_difference())
g_codon_to_missense_codons = dict(_gen_codon_and_missense_codons())


def _aa_pair_to_snp_count(aa_pair):
    snp_count = 0
    aaa, aab = aa_pair
    for ca in g_aa_letter_to_codons[aaa]:
        for cb in g_aa_letter_to_codons[aab]:
            if _hamming_distance(ca, cb) == 1:
                snp_count += 1
    return snp_count
            


def _gen_all_ordered_aa_letter_pairs():
    for a in g_aa_letters:
        for b in g_aa_letters:
            yield (a, b)

g_aa_pair_to_snp_count = dict((pair, _aa_pair_to_snp_count(pair)) for pair in _gen_all_ordered_aa_letter_pairs())

def codon_distribution_to_aa_distribution(codon_distribution):
    """
    @param codon_distribution: maps codons to frequencies
    @return: a mapping from amino acid letters to frequencies
    """
    if set(codon_distribution) != set(g_non_stop_codons):
        raise CodonError('expected a frequency for each codon')
    aa_distribution = {}
    for aa, codons in g_aa_letter_to_codons.items():
        aa_distribution[aa] = sum(codon_distribution[codon] for codon in codons)
    return aa_distribution

def codon_distribution_to_nt_distribution(codon_distribution):
    """
    @param codon_distribution: maps codons to frequencies
    @return: a mapping from amino acid letters to frequencies
    """
    if set(codon_distribution) != set(g_non_stop_codons):
        raise CodonError('expected a frequency for each codon')
    nt_distribution = {'A':0, 'C':0, 'G':0, 'T':0}
    for codon, codon_frequency in codon_distribution.items():
        for nt in codon:
            nt_distribution[nt] += codon_frequency / 3
    return nt_distribution


class Test(unittest.TestCase):
    def test_could_be_missense(self):
        """
        Tests a few amino acid pairs.
        The minimum number of codon nucleotide changes between
        amino acid pairs in the first list is equal to 1.
        This is not true for the second list.
        """
        expected_true = ('AS', 'AP', 'AT', 'AV', 'CW', 'AV', 'AD', 'AE', 'AG')
        expected_false = ('AA', 'AC', 'AY', 'AH', 'AQ', 'AN', 'AK', 'FG', 'FE')
        for pair in expected_true:
            self.assertEqual(g_aa_letter_pair_to_min_codon_difference[tuple(pair)], 1)
        for pair in expected_false:
            self.assertNotEqual(g_aa_letter_pair_to_min_codon_difference[tuple(pair)], 1)


if __name__ == '__main__':
    unittest.main()

