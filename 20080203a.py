"""Given a first order Markov model, generate an endpoint-constrained N-character string.

In the transition matrix, the 27 states are ordered from space to z.
The default transition matrix was estimated from A Tale of Two Cities by Charles Dickens.
"""

import StringIO

import numpy

from SnippetUtil import HandlingError
import PathSampler
import EnglishModel
import Util
import MatrixUtil
import Form

def get_form():
    """
    @return: the body of a form
    """
    # define the default transition matrix
    dictionary_rate_matrix = EnglishModel.get_transition_matrix()
    labels = list(sorted(set(a for a, b in dictionary_rate_matrix)))
    T = numpy.array(MatrixUtil.dict_to_row_major(dictionary_rate_matrix, labels, labels))
    # define the form objects
    form_objects = [
            Form.Matrix('matrix', 'transition matrix', T, MatrixUtil.assert_transition_matrix),
            Form.SingleLine('first', 'first letter', 'a'),
            Form.SingleLine('last', 'last letter', 'b'),
            Form.Integer('count', 'character count', 10, low=1, high=80)]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # define the states using the default transition matrix
    default_transition_matrix = EnglishModel.get_transition_matrix()
    states = list(sorted(set(a for a, b in default_transition_matrix)))
    # read the constraints from the form data
    if fs.first not in states:
        raise HandlingError('invalid first letter')
    if fs.last not in states:
        raise HandlingError('invalid last letter')
    if fs.count == 1 and fs.first != fs.last:
        raise HandlingError('if a single letter is to be simulated then the first letter must be the same as the last letter.')
    # read the transition matrix from the form data
    T = fs.matrix
    if T.shape[0] != len(states):
        raise HandlingError('expected the transition matrix to have %d lines' % len(states))
    matrix = MatrixUtil.row_major_to_dict(T.tolist(), states, states)
    # simulate the path
    path = PathSampler.get_discrete_path_sample(fs.first, fs.last, states, fs.count, matrix)
    # show the simulated path in convenient text form
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, ''.join(path)
