"""Given a transition matrix, get the stationary distribution.
"""

import numpy as np

from SnippetUtil import HandlingError
import SnippetUtil
import EnglishModel
import MatrixUtil
import TransitionMatrix
import Form
import FormOut

#FIXME numpy may not be necessary

def get_form():
    """
    @return: the body of a form
    """
    # define the default transition matrix
    dictionary_rate_matrix = EnglishModel.get_transition_matrix()
    labels = list(sorted(set(a for a, b in dictionary_rate_matrix)))
    T = MatrixUtil.dict_to_row_major(dictionary_rate_matrix, labels, labels)
    T = np.array(T)
    form_objects = [
            Form.Matrix('matrix', 'transition matrix',
                T, MatrixUtil.assert_transition_matrix)]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):
    # get the stationary distribution of the transition matrix
    T = fs.matrix
    try:
        v = TransitionMatrix.get_stationary_distribution(T.tolist())
    except TransitionMatrix.TransitionMatrixError as e
        raise HandlingError(e)
    # return the stationary distribution
    return '\n'.join(str(x) for x in v) + '\n'
