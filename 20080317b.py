"""Analyze a nt mutation distribution and an aa stationary distribution.

Analyze a nucleotide mutation distribution
and an amino acid stationary distribution.
Calculate two output vectors from these input distributions.
The first output vector is the stationary nucleotide distribution.
The second output vector is the centered amino acid energy vector.
"""

from StringIO import StringIO
import math

from SnippetUtil import HandlingError
import SnippetUtil
import Codon
import DirectProtein
import Form
import FormOut
from Codon import g_sorted_nt_letters as nt_letters
from Codon import g_sorted_aa_letters as aa_letters
from Codon import g_sorted_non_stop_codons as codons

def get_form():
    """
    @return: the body of a form
    """
    # define some default strings
    default_nt_string = '\n'.join(nt + ' : 1' for nt in nt_letters)
    default_aa_string = '\n'.join(aa + ' : 1' for aa in aa_letters)
    # define the form objects
    form_objects = [
            Form.MultiLine('nucleotides', 'nucleotide mutation distribution',
                default_nt_string),
            Form.MultiLine('aminoacids', 'amino acid stationary distribution',
                default_aa_string)]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):
    # get the nucleotide distribution
    nt_to_weight = SnippetUtil.get_distribution(fs.nucleotides,
            'nucleotide', nt_letters)
    # get the amino acid distribution
    aa_to_weight = SnippetUtil.get_distribution(fs.aminoacids,
            'amino acid', aa_letters)
    # get results
    mutation_distribution = [nt_to_weight[nt] for nt in nt_letters]
    aa_distribution = [aa_to_weight[aa] for aa in aa_letters]
    pair = DirectProtein.get_nt_distribution_and_aa_energies(
            mutation_distribution, aa_distribution)
    nt_distribution, aa_energies = pair
    # write something
    out = StringIO()
    # write the stationary nucleotide distribution
    print >> out, 'nucleotide stationary distribution:'
    for nt, value in zip(nt_letters, nt_distribution):
        print >> out, '%s : %s' % (nt, value)
    print >> out, ''
    # write the amino acid energies
    print >> out, 'amino acid energies:'
    for aa, value in zip(aa_letters, aa_energies):
        print >> out, '%s : %s' % (aa, value)
    # return the response
    return out.getvalue()
