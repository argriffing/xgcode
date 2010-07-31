"""Estimate F84 parameters from two nt sequences by numerical ML.

Estimate F84 parameters from a pair of nucleotide sequences
by numerically maximizing likelihood.
The F84 evolutionary model is defined in the paper
"Maximum Likelihood Phylogenetic Estimation from DNA Sequences
with Variable Rates over Sites: Approximate Methods"
by Ziheng Yang in J Mol Evol 1994.
"""

from StringIO import StringIO

import scipy.optimize

from SnippetUtil import HandlingError
import Fasta
import F84
import PairLikelihood
import Form
import FormOut
import const

g_sample_fasta_string = const.read('20100730z')

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('fasta', 'nucleotide sequence pair',
                g_sample_fasta_string.strip())]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # get the alignment object
    try:
        alignment = Fasta.Alignment(StringIO(fs.fasta))
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
        msg = 'expected a gapless unambiguous nucleotide alignment'
        raise HandlingError(msg)
    # get the maximum likelihood estimates according to a numeric optimizer.
    f = F84.Objective(alignment.sequences)
    values = list(f.get_initial_parameters())
    result = scipy.optimize.fmin(f, values, ftol=1e-10, disp=0)
    distance, kappa, wC, wG, wT= result
    nt_distribution = F84.parameters_to_distribution((wC, wG, wT))
    A, C, G, T = nt_distribution
    model = F84.create_rate_matrix(kappa, nt_distribution)
    log_likelihood = PairLikelihood.get_log_likelihood(
            distance, alignment.sequences, model)
    # begin the response
    out = StringIO()
    print >> out, 'ML distance:', distance
    print >> out, 'ML kappa:', kappa
    print >> out, 'ML A frequency:', A
    print >> out, 'ML C frequency:', C
    print >> out, 'ML G frequency:', G
    print >> out, 'ML T frequency:', T
    print >> out, 'log likelihood:', log_likelihood
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
