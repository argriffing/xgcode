"""Find the matrix of ML distances between pairs of sequences using JC69.

Find the matrix of maximum likelihood distances
between pairs of aligned sequences using JC69.
The rows and columns of the output distance matrix
are ordered according to the sequence order in the alignment.
The likelihood maximization can be done analytically
because of the simplicity of the model.
"""

import math
from StringIO import StringIO

from SnippetUtil import HandlingError
import Fasta
import MatrixUtil
import JC69
import Form
import FormOut

#FIXME use const data

# HKY simulation parameters:
# transition/transversion ratio 2, C:4, G:4, A:1, T:1, scaled to one
g_fasta = """
>a
TGGGGGCGCACTGTCGGCCGCGTTGCGCAGCACTCCACCGCCTGGCGCTCCCCCGGCGCC
GTGCGGGGGTCGCCGGCGAGCGCGGCTTCCCCAGCGGGCT
>b
GGTGTCCGCTAGGCAGGTGCGTCGCCCCACGTCCGCCCGCGGCGGCCCCGGGACGGCCAC
CCGGGGGAGGGTCGGCGCGACTCGATGCGGTGCCCATGGC
>c
GAATCGCCTCGCGCACCGACGCCGGGCCGGCGTGGCGCCTATCGCGCCCTCGGATGTGTA
CGACGCGCCCGACGCCGTCCCCCTCGCGCGGGGCGCCCTC
>d
CGTCCCTACGCGCGTCCGGGGCGTCGAGCGCGAGCCAGGGCGAATGACTGGCGGATCGCC
GCGGGAGCGGCTCCGTCTAGCCGCACGGGGCGCCCCGTCC
"""

def get_form():
    """
    @return: the body of a form
    """
    # define the form objects
    form_objects = [
            Form.MultiLine('fasta', 'aligned sequences without gaps',
                g_fasta.strip()),
            Form.Float('infinity', 'use this value for estimates of infinity',
                100.0)]
    return form_objects

def get_form_out():
    return FormOut.Matrix()

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # read the alignment
    try:
        alignment = Fasta.Alignment(StringIO(fs.fasta))
    except Fasta.AlignmentError, e:
        raise HandlingError('fasta alignment error: ' + str(e))
    if alignment.get_sequence_count() < 2:
        raise HandlingError('expected at least two sequences')
    # Create the distance matrix,
    # replacing values of None with the representation for infinity.
    row_major_distance_matrix = []
    for row in JC69.get_ML_distance_matrix(alignment.sequences):
        corrected_row = [fs.infinity if x == float('inf') else x for x in row]
        row_major_distance_matrix.append(corrected_row)
    # write the response
    out = StringIO()
    print >> out, MatrixUtil.m_to_string(row_major_distance_matrix)
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
