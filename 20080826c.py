"""Count the number of observed differences between two nucleotide sequences.
"""

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
        msg = 'expected a gapless unambiguous nucleotide alignment'
        raise HandlingError(msg)
    # get the maximum likelihood estimate
    sequence_pair = alignment.sequences
    count = sum(a != b for a, b in zip(*alignment.sequences))
    # return the response
    return 'difference count: %d\n' % count
