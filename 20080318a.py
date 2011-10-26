"""Show the symmetric rate components for a Direct Protein mixture model.
"""

import math
from StringIO import StringIO

from SnippetUtil import HandlingError
import DirectProtein
import RateMatrix
import Form
import FormOut
from Codon import g_sorted_nt_letters as nt_letters
from Codon import g_sorted_aa_letters as aa_letters
from Codon import g_sorted_non_stop_codons as codons

def get_form():
    """
    @return: the body of a form
    """
    default_xml_string = DirectProtein.get_sample_xml_string().strip()
    return [Form.MultiLine('model', 'mixture model', default_xml_string)]

def get_form_out():
    return FormOut.Html()

def get_response_content(fs):
    # deserialize the xml data to create a DirectProteinMixture
    try:
        mixture_model = DirectProtein.deserialize_mixture_model(fs.model)
    except ValueError as e
        raise HandlingError(e)
    # Normalize the mixture model to have an expected rate of one
    # substitution per unit of branch length.
    mixture_model.normalize()
    # begin writing the html file
    out = StringIO()
    # write the html header
    print >> out, '<html>'
    print >> out, '<head>'
    print >> out, '<style type="text/css">td{font-size:x-small;}</style>'
    print >> out, '</head>'
    print >> out, '<body>'
    # write the symmetric components of the rate matrices
    for category_i, matrix_object in enumerate(mixture_model.rate_matrices):
        codon_v = matrix_object.get_stationary_distribution()
        matrix = matrix_object.dictionary_rate_matrix
        symmetric_matrix = {}
        for ca, pa in zip(codons, codon_v):
            for cb, pb in zip(codons, codon_v):
                value = matrix[(ca, cb)] / (math.sqrt(pb) / math.sqrt(pa))
                symmetric_matrix[(ca, cb)] = value
        print >> out, 'the symmetric component of the rate matrix'
        print >> out, 'for category %d:' % (category_i + 1)
        print >> out, '<table>'
        print >> out, RateMatrix.codon_rate_matrix_to_html_string(
                symmetric_matrix)
        print >> out, '</table>'
        print >> out, '<br/><br/>'
    # write the html footer
    print >> out, '</body>'
    print >> out, '</html>'
    # return the response
    return out.getvalue()
