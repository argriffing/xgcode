"""
Given a nt distribution and an aa distribution, get a codon distribution.

Given a nucletoide distribution and an amino acid distribution,
get a codon distribution.
Calculate codon frequencies according to equation (14) of Halpern-Bruno 1998.
"""

from StringIO import StringIO

import numpy as np

from SnippetUtil import HandlingError
import SnippetUtil
import Codon
import Form
import FormOut
from Codon import g_sorted_nt_letters as nt_letters
from Codon import g_sorted_aa_letters as aa_letters

def get_form():
    """
    @return: the body of a form
    """
    # define the default nucleotide and amino acid strings
    default_nt_string = '\n'.join(nt + ' : 1' for nt in nt_letters)
    default_aa_string = '\n'.join(aa + ' : 1' for aa in aa_letters)
    # define the form objects
    form_objects = [
            Form.MultiLine('nucleotides', 'nucleotide distribution weights',
                default_nt_string),
            Form.MultiLine('amino_acids', 'amino acid distribution weights',
                default_aa_string)]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):
    # get the nucleotide distribution
    nt_distribution = SnippetUtil.get_distribution(fs.nucleotides,
            'nucleotide', Codon.g_nt_letters)
    # get the amino acid distribution
    aa_distribution = SnippetUtil.get_distribution(fs.amino_acids,
            'amino acid', Codon.g_aa_letters)
    # Assert that the nucleotide distribution
    # is compatible with the amino acid distribution.
    # According to the Halpern-Bruno assumptions, there should be no codon bias.
    # This means that if a nucleotide has a frequency of zero,
    # then the amino acid coded by each codon containing that nucleotide
    # must also have a frequency of zero.
    msg_a = 'the given amino acid and nucleotide distributions '
    msg_b = 'are incompatible with the assumption of no codon bias'
    err = HandlingError(msg_a + msg_b)
    for aa, codons in Codon.g_aa_letter_to_codons.items():
        for codon in codons:
            for nt in codon:
                if aa_distribution[aa] and not nt_distribution[nt]:
                    raise err
    # get the codon distribution
    codon_to_weight = {}
    for codon in Codon.g_non_stop_codons:
        aa = Codon.g_codon_to_aa_letter[codon]
        sibling_codons = Codon.g_aa_letter_to_codons[aa]
        codon_aa_weight = aa_distribution[aa]
        codon_nt_weight = np.prod([nt_distribution[nt] for nt in codon])
        sibling_nt_weight_sum = 0
        for sibling in sibling_codons:
            product = np.prod([nt_distribution[nt] for nt in sibling])
            sibling_nt_weight_sum += product
        codon_to_weight[codon] = codon_aa_weight * codon_nt_weight
        codon_to_weight[codon] /= sibling_nt_weight_sum
    total_weight = sum(codon_to_weight.values())
    # return the codon distribution
    out = StringIO()
    for codon, weight in sorted(codon_to_weight.items()):
        print >> out, codon, ':', weight / total_weight
    return out.getvalue() + '\n'
