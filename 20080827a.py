"""Find the matrix of ML distances between pairs of sequences using JC69.

Find the matrix of maximum likelihood distances
between pairs of aligned sequences using JC69.
The rows and columns of the output distance matrix
are ordered according to the sequence order in the alignment.
The likelihood maximization can be done analytically
because of the simplicity of the model.
"""

from SnippetUtil import HandlingError
import Fasta
import MatrixUtil
import JC69
import Form
import FormOut
import const

g_fasta = const.read('20100730y')

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

def get_response_content(fs):
    # read the alignment
    try:
        alignment = Fasta.Alignment(fs.fasta.splitlines())
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
    # return the response
    return MatrixUtil.m_to_string(row_major_distance_matrix) + '\n'
