"""Translate a list of codons to amino acids.
"""

from StringIO import StringIO

from SnippetUtil import HandlingError
import Codon
import Form
import FormOut
import Util

def get_form():
    """
    @return: the body of a form
    """
    # define the default list of codons
    default_codons = [
            'TTG', 'CTG', 'TTG', 'CTG', 'CTG',
            'CTG', 'CTC', 'TTG', 'CTG', 'CTG', 'TTG', 'TTG', 'TTG', 'CTG']
    # define the form objects
    form_objects = [
            Form.MultiLine('codons', 'one codon on each line',
                '\n'.join(default_codons))]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # read the list of codons
    codons = Util.get_stripped_lines(StringIO(fs.codons))
    # convert codons to upper case
    codons = [codon.upper() for codon in codons]
    # make sure that each codon is valid
    invalid_codons = set(codon for codon in codons
            if codon not in Codon.g_non_stop_codons)
    if invalid_codons:
        raise HandlingError('invalid codons: ' + ', '.join(invalid_codons))
    # define the response
    out = StringIO()
    for codon in codons:
        aa_letter = Codon.g_codon_to_aa_letter[codon]
        print >> out, aa_letter
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
