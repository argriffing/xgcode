"""Show the symmetric rate components for a Direct Protein mixture model.
"""

import math
import StringIO

from SnippetUtil import HandlingError
import Util
import DirectProtein
import RateMatrix
import Form
from Codon import g_sorted_nt_letters as nt_letters
from Codon import g_sorted_aa_letters as aa_letters
from Codon import g_sorted_non_stop_codons as codons

def get_form():
    """
    @return: the body of a form
    """
    default_xml_string = DirectProtein.get_sample_xml_string().strip()
    return [Form.MultiLine('model', 'mixture model', default_xml_string)]

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # deserialize the xml data to create a DirectProteinMixture
    try:
        mixture_model = DirectProtein.deserialize_mixture_model(fs.model)
    except ValueError, e:
        raise HandlingError(e)
    # Normalize the mixture model to have an expected rate of one
    # substitution per unit of branch length.
    mixture_model.normalize()
    # begin writing the html file
    out = StringIO.StringIO()
    # write the html header
    print >> out, '<html>'
    print >> out, '<head><style type="text/css">td{font-size:x-small;}</style></head>'
    print >> out, '<body>'
    # write the symmetric components of the rate matrices
    for category_index, matrix_object in enumerate(mixture_model.rate_matrices):
        codon_stationary_distribution = matrix_object.get_stationary_distribution()
        matrix = matrix_object.dictionary_rate_matrix
        symmetric_matrix = {}
        for ca, pa in zip(codons, codon_stationary_distribution):
            for cb, pb in zip(codons, codon_stationary_distribution):
                symmetric_matrix[(ca, cb)] = matrix[(ca, cb)] / (math.sqrt(pb) / math.sqrt(pa))
        print >> out, 'the symmetric component of the rate matrix for category %d:' % (category_index + 1)
        print >> out, '<table>'
        print >> out, RateMatrix.codon_rate_matrix_to_html_string(symmetric_matrix)
        print >> out, '</table>'
        print >> out, '<br/><br/>'
    # write the html footer
    print >> out, '</body>'
    print >> out, '</html>'
    # return the response
    response_headers = [('Content-Type', 'text/html')]
    return response_headers, out.getvalue().strip()
