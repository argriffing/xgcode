"""Given a Markov model, generate an endpoint-constrained N-character string.

The input Markov model is first order.
In the transition matrix, the 27 states are ordered from space to z.
The default transition matrix was estimated from
A Tale of Two Cities by Charles Dickens.
"""

from StringIO import StringIO

import numpy as np

from SnippetUtil import HandlingError
import PathSampler
import EnglishModel
import MatrixUtil
import Form
import FormOut

def get_form():
    """
    @return: the body of a form
    """
    # define the default transition matrix
    dictionary_rate_matrix = EnglishModel.get_transition_matrix()
    labels = list(sorted(set(a for a, b in dictionary_rate_matrix)))
    T = MatrixUtil.dict_to_row_major(dictionary_rate_matrix, labels, labels)
    T = np.array(T)
    # define the form objects
    form_objects = [
            Form.Matrix('matrix', 'transition matrix',
                T, MatrixUtil.assert_transition_matrix),
            Form.SingleLine('first', 'first letter', 'a'),
            Form.SingleLine('last', 'last letter', 'b'),
            Form.Integer('count', 'character count', 10, low=1, high=80)]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):
    # define the states using the default transition matrix
    default_transition_matrix = EnglishModel.get_transition_matrix()
    states = list(sorted(set(a for a, b in default_transition_matrix)))
    # read the constraints from the form data
    if fs.first not in states:
        raise HandlingError('invalid first letter')
    if fs.last not in states:
        raise HandlingError('invalid last letter')
    if fs.count == 1 and fs.first != fs.last:
        msg_a = 'if a single letter is to be simulated '
        msg_b = 'then the first letter must be the same as the last letter.'
        raise HandlingError(msg_a + msg_b)
    # read the transition matrix from the form data
    T = fs.matrix
    if T.shape[0] != len(states):
        msg = 'expected the transition matrix to have %d lines' % len(states)
        raise HandlingError(msg)
    matrix = MatrixUtil.row_major_to_dict(T.tolist(), states, states)
    # simulate the path
    path = PathSampler.get_discrete_path_sample(fs.first, fs.last,
            states, fs.count, matrix)
    # show the simulated path in convenient text form
    return ''.join(path) + '\n'
