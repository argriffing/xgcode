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

import math
from StringIO import StringIO

from SnippetUtil import HandlingError
import SnippetUtil
import Fasta
import JC69
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
    form_objects = [
            Form.MultiLine('fasta', 'nucleotide sequence pair',
                g_sample_fasta_string.strip())]
    return form_objects

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
        msg = 'expected a gap-free unambiguous nucleotide alignment'
        raise HandlingError(msg)
    # get the maximum likelihood estimate
    sequence_pair = alignment.sequences
    mle = JC69.get_ML_distance(*sequence_pair)
    # begin the response
    out = StringIO()
    print >> out, 'ML distance estimate:', mle
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
