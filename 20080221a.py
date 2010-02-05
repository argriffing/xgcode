"""Create a Direct Protein codon rate matrix.

Implement the "Direct Protein Model" described by equations (15) and (21)
in the 2007 MBE paper "Population Genetics Without Intraspecific Data"
by Thorne et al.
The kappa parameter is the mutation process transition transversion rate ratio.
"""

from StringIO import StringIO
import math

from SnippetUtil import HandlingError
import SnippetUtil
import Codon
import MatrixUtil
import DirectProtein
import Form
from Codon import g_sorted_nt_letters as nt_ordered
from Codon import g_sorted_aa_letters as aa_ordered
from Codon import g_sorted_non_stop_codons as codons_ordered

def get_form():
    """
    @return: the body of a form
    """
    # define the default nucleotide and amino acid strings
    default_nucleotide_string = '\n'.join(nt + ' : 1' for nt in nt_ordered)
    default_amino_acid_string = '\n'.join(aa + ' : 0' for aa in aa_ordered)
    # define the form objects
    form_objects = [
            Form.MultiLine('nucleotides', 'mutation process nt weights',
                default_nucleotide_string),
            Form.Float('kappa', 'transition transversion rate ratio kappa',
                1, low_exclusive=0),
            Form.MultiLine('aminoacids', 'aa associated unscaled energies',
                default_amino_acid_string),
            Form.RadioGroup('outputtype', 'output options', [
                Form.RadioItem('srm', 'scaled codon rate matrix', True),
                Form.RadioItem('urm', 'unscaled codon rate matrix'),
                Form.RadioItem('cstat', 'codon stationary distribution'),
                Form.RadioItem('astat', 'amino acid stationary distribution'),
                Form.RadioItem('nstat', 'nucleotide stationary distribution'),
                Form.RadioItem('sf', 'rate matrix scaling factor')])]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # get the mutation process nucleotide distribution
    nt_distribution = SnippetUtil.get_distribution(fs.nucleotides,
            'nucleotide', nt_ordered)
    # get the selection process amino acid energies
    aa_to_energy = SnippetUtil.get_dictionary(fs.aminoacids,
            'amino acid', 'energy', aa_ordered)
    # create the direct protein rate matrix object
    nt_distribution_list = [nt_distribution[nt] for nt in nt_ordered]
    aa_energy_list = [aa_to_energy[aa] for aa in aa_ordered]
    rate_matrix_object = DirectProtein.DirectProteinRateMatrix(fs.kappa,
            nt_distribution_list, aa_energy_list)
    # write the response
    out = StringIO()
    if fs.srm:
        # write the scaled rate matrix
        rate_matrix_object.normalize()
        row_major_rate_matrix = rate_matrix_object.get_row_major_rate_matrix()
        print >> out, MatrixUtil.m_to_string(row_major_rate_matrix)
    elif fs.urm:
        # write the unscaled rate matrix
        row_major_rate_matrix = rate_matrix_object.get_row_major_rate_matrix()
        print >> out, MatrixUtil.m_to_string(row_major_rate_matrix)
    elif fs.cstat:
        # write the codon stationary distribution
        codon_distribution = rate_matrix_object.get_codon_distribution()
        for codon in codons_ordered:
            print >> out, codon, ':', codon_distribution[codon]
    elif fs.astat:
        # write the amino acid stationary distribution
        aa_distribution = rate_matrix_object.get_aa_distribution()
        for aa in aa_ordered:
            print >> out, aa, ':', aa_distribution[aa]
    elif fs.nstat:
        # write the nucleotide stationary distribution
        nt_distribution = rate_matrix_object.get_nt_distribution()
        for nt in nt_ordered:
            print >> out, nt, ':', nt_distribution[nt]
    elif fs.sf:
        # write the rate matrix scaling factor
        print >> out, rate_matrix_object.get_expected_rate()
    # return the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
