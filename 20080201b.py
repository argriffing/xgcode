"""Create a GY94 codon rate matrix from codon weights and two other parameters.

The rows and columns in the output will be ordered alphabetically by codon regardless of the input order.
Kappa controls the transition / transversion rate.
Omega controls the non-synonymous / synonymous rate.
"""

from StringIO import StringIO

from SnippetUtil import HandlingError
import SnippetUtil
import CodonFrequency
import Codon
import RateMatrix
import Form

def get_form():
    """
    @return: the body of a form
    """
    # define the default codon weight lines
    default_weight_lines = []
    cf = CodonFrequency.codon_frequency_b
    for codon in Codon.g_sorted_non_stop_codons:
        weight = cf.codon_to_non_stop_proportion(codon)
        line = codon + ' : ' + str(weight)
        default_weight_lines.append(line)
    # define the form objects
    form_objects = [
            Form.MultiLine('weights', 'codon weights', '\n'.join(default_weight_lines)),
            Form.Float('kappa', 'kappa', 2),
            Form.Float('omega', 'omega', .01)]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # get the codon distribution
    codons = Codon.g_sorted_non_stop_codons
    distribution = SnippetUtil.get_distribution(fs.weights, 'codon', codons)
    # get the rate matrix defined by the weights and kappa and omega
    r = RateMatrix.get_gy94_rate_matrix(distribution, fs.kappa, fs.omega)
    # show the rate matrix in convenient text form
    out = StringIO()
    for ca in codons:
        print >> out, '\t'.join(str(r[(ca, cb)]) for cb in codons)
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
