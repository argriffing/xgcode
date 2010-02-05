"""Given a rate matrix, get the stationary distribution.
"""

from StringIO import StringIO

import numpy

from SnippetUtil import HandlingError
import Util
import MatrixUtil
import RateMatrix
import Form

def get_form():
    """
    @return: the body of a form
    """
    # define the default rate matrix
    dictionary_rate_matrix = RateMatrix.get_sample_codon_rate_matrix()
    labels = list(sorted(set(a for a, b in dictionary_rate_matrix)))
    R = numpy.array(MatrixUtil.dict_to_row_major(dictionary_rate_matrix, labels, labels))
    return [Form.Matrix('matrix', 'rate matrix', R, MatrixUtil.assert_rate_matrix)]

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # read the matrix from the form data
    R = fs.matrix
    # get the stationary distribution of the rate matrix
    try:
        stationary_distribution = RateMatrix.get_stationary_distribution(R.tolist())
    except RateMatrix.RateMatrixError, e:
        raise HandlingError('error calculating the stationary distribution: ' + str(e))
    # get the stationary distribution string
    stationary_distribution_string = '\n'.join(str(x) for x in stationary_distribution)
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, stationary_distribution_string
