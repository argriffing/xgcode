"""Analyze a single nucleotide site of a pileup.

Parameters x, y, and z are branch lengths which define the two three-taxon
trees whose Jukes-Cantor nucleotide distribution at the tips
define the {RR, RA, AA, AB} zygosity distributions for the
recent vs. ancient mrca states.
The low, medium, and high parameters are expectations of three geometrically
distributed mixture components of a garbage state.
The seqerror parameter is the probability of sequencing randomization;
this the probability that the sequencing machine spits out a random
nucleotide instead of the correct nucleotide.
The nomcoverage parameter defines the nominal coverage of the pileup.
The kmulticoverages parameter defines the number of
nominal coverage multiples which might result from duplications.
"""

from StringIO import StringIO

from SnippetUtil import HandlingError
import Form
import DGRP

def get_form():
    """
    @return: the body of a form
    """
    default_params = [
            ('x', '0.1'),
            ('y', '0.01'),
            ('z', '0.0001'),
            ('low', '2'),
            ('med', '20'),
            ('high', '1000'),
            ('seqerr', '.01'),
            ('nomcoverage', '20'),
            ('kmulticoverages', '4')]
    default_counts = [
            ('A', '20'),
            ('C', '20'),
            ('G', '3'),
            ('T', '1')]
    form_objects = [
            Form.MultiLine('param_field', 'parameters',
                '\n'.join('\t'.join(p) for p in default_params)),
            Form.MultiLine('count_field', 'nucleotide read counts',
                '\n'.join('\t'.join(p) for p in default_counts)),
            Form.SingleLine('ref', 'reference nucleotide', 'A')]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    out = StringIO()
    return [('Content-Type', 'text/plain')], out.getvalue().strip()
