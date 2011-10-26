"""Summarize the properties of a Direct Protein mixture model.
"""

from StringIO import StringIO

from SnippetUtil import HandlingError
import Codon
import DirectProtein
import Form
import FormOut

def get_form():
    """
    @return: the body of a form
    """
    # define the default xml string
    default_xml_string = DirectProtein.get_sample_xml_string()
    # define the form objects
    return [Form.MultiLine('model', 'mixture model', default_xml_string)]

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):
    # deserialize the xml data to create a DirectProteinMixture
    try:
        mixture_model = DirectProtein.deserialize_mixture_model(fs.model)
    except ValueError as e:
        raise HandlingError(e)
    expected_rate = mixture_model.get_expected_rate()
    codon_distribution = mixture_model.get_codon_stationary_distribution()
    aa_distribution = mixture_model.get_aa_stationary_distribution()
    nt_distribution = mixture_model.get_nt_stationary_distribution()
    ordered_codons = list(sorted(Codon.g_non_stop_codons))
    # show the summary
    out = StringIO()
    print >> out, 'expected codon substitution rate:'
    print >> out, expected_rate
    print >> out, ''
    print >> out, 'nucleotide distribution:'
    for nt, proportion in zip(Codon.g_nt_letters, nt_distribution):
        print >> out, '%s : %s' % (nt, proportion)
    print >> out, ''
    print >> out, 'amino acid distribution:'
    for aa, proportion in zip(Codon.g_aa_letters, aa_distribution):
        print >> out, '%s : %s' % (aa, proportion)
    print >> out, ''
    print >> out, 'codon distribution:'
    for codon, proportion in zip(ordered_codons, codon_distribution):
        print >> out, '%s : %s' % (codon, proportion)
    return out.getvalue()
