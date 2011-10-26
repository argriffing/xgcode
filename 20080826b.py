"""Estimate the distance between two sequences using an ML JC69 estimator.

Estimate the distance between two sequences using
an explicit ML JC69 estimator.
JC69 is the simple continuous time Markov model
described by Jukes and Cantor in 1969.
The distance between two sequences is the expected number
of nucleotide changes per position.
For the JC69 model, the maximum likelihood estimator
for the distance can be written explicitly.
"""

from SnippetUtil import HandlingError
import SnippetUtil
import Fasta
import JC69
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

def get_response_content(fs):
    # get the alignment object
    try:
        alignment = Fasta.Alignment(fs.fasta.splitlines())
    except Fasta.AlignmentError as e
        raise HandlingError('alignment error: ' + str(e))
    # assert that the alignment is of exactly two sequences
    if len(alignment.sequences) != 2:
        raise HandlingError('expected a pair of sequences')
    # assert that the alignment is a gapless unambiguous nucleotide alignment
    old_column_count = alignment.get_column_count()
    try:
        alignment.force_nucleotide()
    except Fasta.AlignmentError as e
        raise HandlingError('nucleotide alignment error: ' + str(e))
    new_column_count = alignment.get_column_count()
    if old_column_count != new_column_count:
        msg = 'expected a gap-free unambiguous nucleotide alignment'
        raise HandlingError(msg)
    # get the maximum likelihood estimate
    sequence_pair = alignment.sequences
    mle = JC69.get_ML_distance(*sequence_pair)
    # return the response
    return 'ML distance estimate: %f\n' % mle
