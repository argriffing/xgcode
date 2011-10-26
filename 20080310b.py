"""Summarize the properties of a Direct RNA mixture model.

See "Population Genetics Without Intraspecific Data" by Thorne et al.
for more information about the Direct RNA model.
"""

from StringIO import StringIO

from SnippetUtil import HandlingError
import Codon
import DirectRna
import Form
import FormOut

def get_form():
    """
    @return: the body of a form
    """
    # define the default xml string
    default_xml_string = DirectRna.get_sample_xml_string()
    # define the form objects
    return [Form.MultiLine('model', 'mixture model', default_xml_string)]

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):
    # deserialize the xml data to create a DirectRnaMixture
    try:
        mixture_model = DirectRna.deserialize_mixture_model(fs.model)
    except ValueError as e:
        raise HandlingError(e)
    expected_rate = mixture_model.get_expected_rate()
    nt_distribution = mixture_model.get_nt_stationary_distribution()
    # write the summary
    out = StringIO()
    print >> out, 'expected substitution rate:'
    print >> out, expected_rate
    print >> out, ''
    print >> out, 'nucleotide distribution:'
    for nt, proportion in zip(Codon.g_nt_letters, nt_distribution):
        print >> out, '%s : %s' % (nt, proportion)
    # return the summary
    return out.getvalue()
