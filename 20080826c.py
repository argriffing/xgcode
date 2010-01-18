"""Count the number of observed differences between two nucleotide sequences.
"""

import math
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
    # get the maximum likelihood estimate
    sequence_pair = alignment.sequences
    count = sum(1 for a, b in zip(*alignment.sequences) if a != b)
    # begin the response
    out = StringIO.StringIO()
    print >> out, 'difference count:', count
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
