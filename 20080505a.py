"""Create an HKY-like nucleotide rate matrix with a WAG-like f parameter.

The rows and columns in the output will be ordered alphabetically
by nucleotide (A, C, G, T) regardless of the input order.
For an HKY85 rate matrix set f to 1.
For a Jukes-Cantor rate matrix set each nucleotide weight to 1,
and set f and &kappa; to 1.
"""

from StringIO import StringIO

from SnippetUtil import HandlingError
import SnippetUtil
import Newick
import RateMatrix
import MatrixUtil
import Form
import FormOut

def get_form():
    """
    @return: the body of a form
    """
    # define the default nucleotide weight lines
    default_nt_lines = [
            'A : 1',
            'C : 1',
            'G : 1',
            'T : 1']
    # define the form objects
    form_objects = [
            Form.MultiLine('weights', 'nucleotide weights',
                '\n'.join(default_nt_lines)),
            Form.Float('kappa', 'kappa', 2, low_inclusive=0),
            Form.Float('f', 'f', 0.5, low_inclusive=0, high_inclusive=1),
            Form.RadioGroup('format', 'rate matrix scaling options', [
                Form.RadioItem('scaled', 'scaled to a rate of one', True),
                Form.RadioItem('unscaled', 'unscaled')])]
    return form_objects

def get_form_out():
    return FormOut.NucleotideRateMatrix()

def get_response_content(fs):
    # get the nucleotide distribution
    distribution = SnippetUtil.get_distribution(
            fs.weights, 'nucleotide', list('ACGT'))
    # get the rate matrix defined by the nucleotide distribution and kappa
    rate_matrix_object = create_rate_matrix(distribution, fs.kappa, fs.f)
    if fs.scaled:
        rate_matrix_object.normalize()
    rate_matrix = rate_matrix_object.get_dictionary_rate_matrix()
    # show the rate matrix in convenient text form
    out = StringIO()
    for nta in 'ACGT':
        print >> out, '\t'.join(str(rate_matrix[(nta, ntb)]) for ntb in 'ACGT')
    return out.getvalue()

def create_rate_matrix(distribution, kappa, f):
    """
    The parameter f does not affect the stationary distribution.
    @param distribution: a dictionary mapping a nucleotide to its frequency
    @param kappa: the transition / transversion substitution rate ratio
    @param f: a WAG-like parameter between zero and one
    @return: a nucleotide rate matrix object
    """
    assert len(distribution) == 4
    assert set(distribution) == set('ACGT')
    assert abs(sum(distribution.values()) - 1.0) < .0000001
    # Create the off-diagonal elements of the unscaled rate matrix.
    rate_matrix = {}
    for na, pa in distribution.items():
        for nb, pb in distribution.items():
            if na != nb:
                if f == 1:
                    rate = pb
                else:
                    rate = (pb**f) / (pa**(1-f))
                if na+nb in ('AG', 'GA', 'CT', 'TC'):
                    rate *= kappa
                rate_matrix[(na, nb)] = rate
    # Create the diagonal elements 
    # such that each row in the rate matrix sums to zero.
    for na in distribution:
        rate = sum(rate_matrix[(na, nb)] for nb in distribution if nb != na)
        rate_matrix[(na, na)] = -rate
    # Convert the dictionary rate matrix to a row major rate matrix
    ordered_states = list('ACGT')
    row_major_rate_matrix = MatrixUtil.dict_to_row_major(
            rate_matrix, ordered_states, ordered_states)
    rate_matrix_object = RateMatrix.RateMatrix(
            row_major_rate_matrix, ordered_states)
    return rate_matrix_object
