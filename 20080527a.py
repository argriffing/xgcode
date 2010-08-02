"""Calculate the expected time spent in each state on a path.

Given a rate matrix,
calculate the expected amount of time spent in each state
given known initial and final states
and a known amount of time between the states.
The rate matrix is expected to be time reversible.
"""

from StringIO import StringIO

import numpy as np

from SnippetUtil import HandlingError
import SnippetUtil
import MatrixUtil
import RateMatrix
import Form
import FormOut

def get_form():
    """
    @return: the body of a form
    """
    # define the default matrix
    R = np.array([
        [-3, 1, 1, 1],
        [2, -4, 1, 1],
        [2, 1, -4, 1],
        [2, 1, 1, -4]])
    # define the form objects
    form_objects = [
            Form.Matrix('matrix', 'rate matrix',
                R, MatrixUtil.assert_rate_matrix),
            Form.Integer('initial', 'initial state', 0, low=0),
            Form.Integer('final', 'final state', 1, low=0),
            Form.Float('time', 'time between observations',
                2.0, low_exclusive=0)]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):
    # read the matrix from the form data
    R = fs.matrix
    matrix_size = len(R)
    # enforce constraints on the initial and final states
    if fs.initial >= matrix_size:
        msg = 'the initial state must be less than the size of the matrix'
        raise HandlingError(msg)
    if fs.final >= matrix_size:
        msg = 'the final state must be less than the size of the matrix'
        raise HandlingError(msg)
    # convert the row major rate matrix to a rate matrix object
    arbitrary_states = list(str(x) for x in range(matrix_size))
    rate_matrix_object = RateMatrix.FastRateMatrix(
            R.tolist(), arbitrary_states)
    # get the expected time spent in each state
    expected_wait_times = rate_matrix_object.get_expected_times(
            fs.initial, fs.final, fs.time)
    # create the response text
    out = StringIO()
    for state_index, expected_wait_time in enumerate(expected_wait_times):
        print >> out, state_index, ':', expected_wait_time / fs.time
    # return the response
    return out.getvalue()
