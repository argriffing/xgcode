#!/usr/bin/env python

from StringIO import StringIO
import re

import Codon
import Util

class CodonFrequency:

    def codon_to_count(self, codon):
        """
        @return: the number of times the codon was observed
        """
        raise NotImplementedError

    def codon_to_non_stop_proportion(self, codon):
        """
        @return: the proportion of observations of the codon among all non stop codons
        """
        raise NotImplementedError


# These codon counts are from the Kazusa DNA Research Institute.
# The counts were taken from Homo sapiens coding sequences.
# The counts were made available here:
#     - http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9606
# This was the upstream data source:
#     - NCBI-GenBank Flat File Release 160.0 [June 15 2007]. 

_raw_human_codon_count_string = """
UUU 17.6(714298)  UCU 15.2(618711)  UAU 12.2(495699)  UGU 10.6(430311)
UUC 20.3(824692)  UCC 17.7(718892)  UAC 15.3(622407)  UGC 12.6(513028)
UUA  7.7(311881)  UCA 12.2(496448)  UAA  1.0( 40285)  UGA  1.6( 63237)
UUG 12.9(525688)  UCG  4.4(179419)  UAG  0.8( 32109)  UGG 13.2(535595)

CUU 13.2(536515)  CCU 17.5(713233)  CAU 10.9(441711)  CGU  4.5(184609)
CUC 19.6(796638)  CCC 19.8(804620)  CAC 15.1(613713)  CGC 10.4(423516)
CUA  7.2(290751)  CCA 16.9(688038)  CAA 12.3(501911)  CGA  6.2(250760)
CUG 39.6(1611801)  CCG  6.9(281570)  CAG 34.2(1391973)  CGG 11.4(464485)

AUU 16.0(650473)  ACU 13.1(533609)  AAU 17.0(689701)  AGU 12.1(493429)
AUC 20.8(846466)  ACC 18.9(768147)  AAC 19.1(776603)  AGC 19.5(791383)
AUA  7.5(304565)  ACA 15.1(614523)  AAA 24.4(993621)  AGA 12.2(494682)
AUG 22.0(896005)  ACG  6.1(246105)  AAG 31.9(1295568)  AGG 12.0(486463)

GUU 11.0(448607)  GCU 18.4(750096)  GAU 21.8(885429)  GGU 10.8(437126)
GUC 14.5(588138)  GCC 27.7(1127679)  GAC 25.1(1020595)  GGC 22.2(903565)
GUA  7.1(287712)  GCA 15.8(643471)  GAA 29.0(1177632)  GGA 16.5(669873)
GUG 28.1(1143534)  GCG  7.4(299495)  GAG 39.6(1609975)  GGG 16.5(669768)
"""

class CodonFrequencyA(CodonFrequency):
    """
    Use human codon counts derived from GenBank and compiled by Kazusa Institute.
    """
    def __init__(self):
        self._codon_to_count = {}
        self._init_codon_to_count()
        self._non_stop_codon_count = sum(self._codon_to_count[codon] for codon in Codon.g_non_stop_codons)

    def _init_codon_to_count(self):
        # parse the raw string
        codon_pattern = r'([ACGU][ACGU][ACGU])(.*)\((.*)\)'
        line_pattern = '.*'.join([codon_pattern]*4)
        token_lists = []
        for line in StringIO(_raw_human_codon_count_string):
            line = line.strip()
            if line:
                m = re.search(line_pattern, line)
                token_list = [x.strip() for x in m.groups()]
                assert len(token_list) == 12
                token_lists.append(token_list)
        assert len(token_lists) == 16
        # write the dictionary using the tokens
        for token_list in token_lists:
            for token in token_list:
                assert type(token) == str
            for codon, per_thousand, count in Util.chopped(token_list, 3):
                # validate the codon
                assert len(codon) == 3
                assert set(codon) <= set('ACGU')
                # use the dna sense codon instead of the rna codon
                codon = codon.replace('U', 'T')
                # validate the count
                assert set(count) <= set('0123456789')
                count = int(count)
                assert count > 0
                self._codon_to_count[codon] = count
        assert len(self._codon_to_count) == 64

    def codon_to_count(self, codon):
        """
        @return: the number of times the codon was observed
        """
        return self._codon_to_count[codon]

    def codon_to_non_stop_proportion(self, codon):
        """
        @return: the proportion of observations of the codon among all non stop codons
        """
        return self._codon_to_count[codon] / float(self._non_stop_codon_count)


# These codon proportions are from the supplementary data
# of a manuscript by Stone and Hobolth.

_raw_stone_codon_proportion_string = """
cdn aa  frq
AAA Lys 0.01850976
AAC Asn 0.02562776
AAG Lys 0.03957337
AAT Asn 0.02201664
ACA Thr 0.01078495
ACC Thr 0.02062059
ACG Thr 0.01387488
ACT Thr 0.0098617
AGA Arg 0.00556931
AGC Ser 0.01993187
AGG Arg 0.00631759
AGT Ser 0.0120507
ATA Ile 0.01025631
ATC Ile 0.02172254
ATG Met 0.02409396
ATT Ile 0.01745249
CAA Gln 0.01639521
CAC His 0.01505873
CAG Gln 0.03532565
CAT His 0.01104927
CCA Pro 0.01356961
CCC Pro 0.01725146
CCG Pro 0.01479069
CCT Pro 0.00702865
CGA Arg 0.00910969
CGC Arg 0.01710254
CGG Arg 0.00835024
CGT Arg 0.00920276
CTA Leu 0.00911714
CTC Leu 0.01369618
CTG Leu 0.03898889
CTT Leu 0.01005156
GAA Glu 0.0222549
GAC Asp 0.02397483
GAG Glu 0.04255533
GAT Asp 0.02920537
GCA Ala 0.01300374
GCC Ala 0.0311375
GCG Ala 0.01275804
GCT Ala 0.01474601
GGA Gly 0.01655902
GGC Gly 0.02404929
GGG Gly 0.00421421
GGT Gly 0.01264635
GTA Val 0.00657074
GTC Val 0.01278782
GTG Val 0.02835657
GTT Val 0.01173427
TAA Stp 0
TAG Stp 0
TAC Tyr 0.01944419
TAT Tyr 0.0126054
TCA Ser 0.00764663
TCC Ser 0.01899373
TCG Ser 0.01591125
TCT Ser 0.0068909
TGA Stp 0
TGC Cys 0.01253839
TGG Trp 0.01048341
TGT Cys 0.00522309
TTA Leu 0.00448597
TTC Phe 0.02221767
TTG Leu 0.01783966
TTT Phe 0.01481302
"""

class CodonFrequencyB(CodonFrequency):

    def __init__(self):
        self._codon_to_non_stop_proportion = {}
        self._init_codon_to_non_stop_proportion()

    def _init_codon_to_non_stop_proportion(self):
        pattern = r'([ACGT][ACGT][ACGT]) ([A-Z][a-z][a-z]) (.*)'
        for line in StringIO(_raw_stone_codon_proportion_string):
            line = line.strip()
            if not line:
                continue
            m = re.search(pattern, line)
            if m:
                codon, aa_name, proportion = m.groups()
                proportion = float(proportion)
                self._codon_to_non_stop_proportion[codon] = proportion

    def codon_to_non_stop_proportion(self, codon):
        """
        @return: the proportion of observations of the codon among all non stop codons
        """
        return self._codon_to_non_stop_proportion[codon]

# codon frequency object for human codons from genbank
codon_frequency_a = CodonFrequencyA()

# codon frequency object for codons from the manuscript by Stone and Hobolth
codon_frequency_b = CodonFrequencyB()

def main():
    print '...'

if __name__ == '__main__':
    main()

