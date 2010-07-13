"""Use the pruning algorithm to get JC log likelihood of an nt alignment.

Use the pruning algorithm to calculate the logarithm
of the Jukes-Cantor likelihood of a nucleotide alignment.
"""

from StringIO import StringIO

from SnippetUtil import HandlingError
import Newick
import Fasta
import PhyLikelihood
import RateMatrix
import MatrixUtil
import Form
import FormOut

def get_form():
    """
    @return: the body of a form
    """
    # define the tree and alignment strings
    tree_string = Newick.brown_example_tree
    tree = Newick.parse(tree_string, Newick.NewickTree)
    formatted_tree_string = Newick.get_narrow_newick_string(tree, 60)
    alignment_string = Fasta.brown_example_alignment.strip()
    # define the form objects
    form_objects = [
            Form.MultiLine('tree', 'newick tree', formatted_tree_string),
            Form.MultiLine('fasta', 'fasta alignment', alignment_string)]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # get the tree
    tree = Newick.parse(fs.tree, Newick.NewickTree)
    tree.assert_valid()
    # get the alignment
    try:
        alignment = Fasta.Alignment(StringIO(fs.fasta))
        alignment.force_nucleotide()
    except Fasta.AlignmentError, e:
        raise HandlingError(e)
    # get the log likelihood
    dictionary_rate_matrix = RateMatrix.get_jukes_cantor_rate_matrix()
    ordered_states = list('ACGT')
    row_major_rate_matrix = MatrixUtil.dict_to_row_major(
            dictionary_rate_matrix, ordered_states, ordered_states)
    rate_matrix_object = RateMatrix.RateMatrix(
            row_major_rate_matrix, ordered_states)
    log_likelihood = PhyLikelihood.get_log_likelihood(
            tree, alignment, rate_matrix_object)
    # write the response
    out = StringIO()
    print >> out, log_likelihood
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
