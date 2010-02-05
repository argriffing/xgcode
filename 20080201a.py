"""Create an HKY85 nt rate matrix from nt weights and a kappa parameter.

Here nt means nucleotide.
The rows and columns in the output will be ordered alphabetically by nucleotide (A, C, G, T)
regardless of the input order.
For the Jukes-Cantor rate matrix set each nucleotide weight to 1, and set kappa to 1.
"""

from StringIO import StringIO

from SnippetUtil import HandlingError
import SnippetUtil
import Newick
import RateMatrix
import MatrixUtil
import Form

def get_form():
    """
    @return: the body of a form
    """
    # define the default nucleotide weights
    default_nt_weight_lines = [
            'A : 1',
            'C : 1',
            'G : 1',
            'T : 1']
    # define the form objects
    form_objects = [
            Form.MultiLine('weights', 'nucleotide weights',
                '\n'.join(default_nt_weight_lines)),
            Form.Float('kappa', 'kappa', 2),
            Form.RadioGroup('format', 'rate matrix scaling', [
                Form.RadioItem('scaled', 'scaled to a rate of one', True),
                Form.RadioItem('unscaled', 'unscaled')])]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # get the nucleotide distribution
    d = SnippetUtil.get_distribution(fs.weights, 'nucleotide', list('ACGT'))
    # get the rate matrix defined by the nucleotide distribution and kappa
    rate_object = RateMatrix.get_unscaled_hky85_rate_matrix(d, fs.kappa)
    if fs.scaled:
        rate_object.normalize()
    rate_matrix = rate_object.get_dictionary_rate_matrix()
    # show the rate matrix in convenient text form
    out = StringIO()
    for nta in 'ACGT':
        print >> out, '\t'.join(str(rate_matrix[(nta, ntb)]) for ntb in 'ACGT')
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
