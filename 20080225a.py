"""Given a transition matrix, get the stationary distribution.
"""

from StringIO import StringIO

import numpy

from SnippetUtil import HandlingError
import SnippetUtil
import EnglishModel
import MatrixUtil
import TransitionMatrix
import Form

def get_form():
    """
    @return: the body of a form
    """
    # define the default transition matrix
    dictionary_rate_matrix = EnglishModel.get_transition_matrix()
    labels = list(sorted(set(a for a, b in dictionary_rate_matrix)))
    T = numpy.array(MatrixUtil.dict_to_row_major(dictionary_rate_matrix, labels, labels))
    return [Form.Matrix('matrix', 'transition matrix', T, MatrixUtil.assert_transition_matrix)]

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # get the stationary distribution of the transition matrix
    T = fs.matrix
    try:
        stationary_distribution = TransitionMatrix.get_stationary_distribution(T.tolist())
    except TransitionMatrix.TransitionMatrixError, e:
        raise HandlingError(e)
    # get the stationary distribution string
    stationary_distribution_string = '\n'.join(str(x) for x in stationary_distribution)
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, stationary_distribution_string
