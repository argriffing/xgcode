"""Find the ML evolutionary distance between two sequences given a rate matrix.

Find the maximum likelihood evolutionary distance
between a pair of sequences given a rate matrix.
"""

from StringIO import StringIO

import scipy.optimize
import numpy as np

from SnippetUtil import HandlingError
import Fasta
import MatrixUtil
import RateMatrix
import Util
import PairLikelihood
import Form
import FormOut

#FIXME use const data

# HKY simulation parameters:
# transition/transversion ratio 2, C:4, G:4, A:1, T:1, scaled to one
g_fasta = """
>hello
GGCACCGCGGGCGGCCCCGGAGCAGACCGGGTACACGCGTCCGTGCTGCCTCCCGCCCGT
AGGTGGGGGGGCCGCCGTGCCGATGACCACCGGCAGCGATCCCACCGAGAACCGGTCGCC
CCTTGTGCCGATGCCACAGCCAGACAGCCTACCCCGGCTGCTGCGCCGCTCGCGACCGCT
GACGTGAAGGCCGCCGTTTG
>world
CCCCGGGCGACCCCAACCACCGCCGGCCCACCGGCCGCCGGCCTGGCTGGGGGGGAGACG
ATGGGGCTGGTACTCGCCGGGCTAGACCAGCCGGCGTGCCCCCTGCCACAGTCGGCTGTC
CTTGGTCAAGCCAGGCCGAGTCACCAGCCCGAGCCGGGGGGCGCCGCCCGCCCCCCCGCA
CGGCCGGCGGCCAGGGCGGG
"""

def get_form():
    """
    @return: the body of a form
    """
    # define the default matrix
    # HKY simulation parameters are
    # transition/transversion ratio 2, C:4, G:4, A:1, T:1
    R = np.array([
            [-1.3, 0.4, 0.8, 0.1],
            [0.1, -0.7, 0.4, 0.2],
            [0.2, 0.4, -0.7, 0.1],
            [0.1, 0.8, 0.4, -1.3]])
    # define the form objects
    form_objects = [
            Form.MultiLine('fasta', 'a gap-free aligned nucleotide pair',
                g_fasta.strip()),
            Form.Matrix('matrix', 'rate matrix',
                R, MatrixUtil.assert_rate_matrix),
            Form.MultiLine('states', 'ordered_states',
                '\n'.join('ACGT'))]
    return form_objects

def get_form_out():
    return FormOut.Report()


class Objective:

    def __init__(self, sequence_pair, rate_matrix):
        """
        @param sequence_pair: a pair of sequences
        @param rate_matrix: a rate matrix object
        """
        self.rate_matrix = rate_matrix
        self.sequence_pair = sequence_pair

    def __call__(self, branch_length):
        """
        This will be called by a one dimensional minimizer.
        @param branch_length: the distance between the two aligned sequences
        @return: the negative log likelihood
        """
        if branch_length < 0:
            return float('inf')
        log_likelihood = PairLikelihood.get_log_likelihood(
                branch_length, self.sequence_pair, self.rate_matrix)
        if log_likelihood is None:
            return float('inf')
        return -log_likelihood


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
    if alignment.get_sequence_count() != 2:
        raise HandlingError('expected a sequence pair')
    # read the rate matrix
    R = fs.matrix
    # read the ordered states
    ordered_states = Util.get_stripped_lines(StringIO(fs.states))
    if len(ordered_states) != len(R):
        msg_a = 'the number of ordered states must be the same '
        msg_b = 'as the number of rows in the rate matrix'
        raise HandlingError(msg_a + msg_b)
    if len(set(ordered_states)) != len(ordered_states):
        raise HandlingError('the ordered states must be unique')
    # create the rate matrix object using the ordered states
    rate_matrix_object = RateMatrix.RateMatrix(R.tolist(), ordered_states) 
    # create the objective function
    objective = Objective(alignment.sequences, rate_matrix_object)
    # Use golden section search to find the mle distance.
    # The bracket is just a suggestion.
    bracket = (0.51, 2.01)
    mle_distance = scipy.optimize.golden(objective, brack=bracket)
    # write the response
    out = StringIO()
    print >> out, 'maximum likelihood distance:', mle_distance
    #distances = (mle_distance, 0.2, 2.0, 20.0)
    #for distance in distances:
        #print >> out, 'f(%s): %s' % (distance, objective(distance))
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
