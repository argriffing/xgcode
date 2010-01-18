"""Estimate F84 parameters from a pair of nucleotide sequences using a closed form estimator.

The F84 evolutionary model and the closed form estimator are defined in the paper
"Maximum Likelihood Phylogenetic Estimation from DNA Sequences with Variable Rates over Sites: Approximate Methods"
by Ziheng Yang in J Mol Evol 1994.
"""

import StringIO

from SnippetUtil import HandlingError
import Fasta
import F84
import PairLikelihood
import Form

g_sample_fasta_string = """
>sequence_a
AAAACCCCGGGGTTAA
>sequence_b
GAAACCTCGGCGTAAA
"""

def get_form():
    """
    @return: the body of a form
    """
    return [Form.MultiLine('fasta', 'nucleotide sequence pair', g_sample_fasta_string.strip())]

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # get the alignment object
    try:
        alignment = Fasta.Alignment(StringIO.StringIO(fs.fasta))
    except Fasta.AlignmentError, e:
        raise HandlingError('alignment error: ' + str(e))
    # assert that the alignment is of exactly two sequences
    if len(alignment.sequences) != 2:
        raise HandlingError('expected a pair of sequences')
    # assert that the alignment is a gapless unambiguous nucleotide alignment
    old_column_count = alignment.get_column_count()
    try:
        alignment.force_nucleotide()
    except Fasta.AlignmentError, e:
        raise HandlingError('nucleotide alignment error: ' + str(e))
    new_column_count = alignment.get_column_count()
    if old_column_count != new_column_count:
        raise HandlingError('expected a gapless unambiguous nucleotide alignment')
    # get the maximum likelihood estimates
    sequence_pair = alignment.sequences
    distance, kappa, A, C, G, T = F84.get_closed_form_estimates(sequence_pair)
    # get the log likelihood
    nt_distribution = (A, C, G, T)
    rate_matrix_object = F84.create_rate_matrix(kappa, nt_distribution)
    log_likelihood = PairLikelihood.get_log_likelihood(distance, alignment.sequences, rate_matrix_object)
    # begin the response
    out = StringIO.StringIO()
    print >> out, 'distance:', distance
    print >> out, 'kappa:', kappa
    print >> out, 'A frequency:', A
    print >> out, 'C frequency:', C
    print >> out, 'G frequency:', G
    print >> out, 'T frequency:', T
    print >> out, 'log likelihood:', log_likelihood
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
