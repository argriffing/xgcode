"""Given a rate matrix, get the expected number of transitions per unit time.
"""

from StringIO import StringIO

import numpy as np

from SnippetUtil import HandlingError
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
    R = MatrixUtil.dict_to_row_major(dictionary_rate_matrix, labels, labels)
    R = np.array(R)
    form_objects = [
            Form.Matrix('matrix', 'rate matrix',
                R, MatrixUtil.assert_rate_matrix)]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # read the matrix from the form data
    R = fs.matrix
    # get the expected rate
    states = range(len(R))
    try:
        rate_matrix_object = RateMatrix.RateMatrix(R.tolist(), states)
        expected_rate = rate_matrix_object.get_expected_rate()
    except RateMatrix.RateMatrixError, e:
        raise HandlingError('error calculating the expected rate: ' + str(e))
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, str(expected_rate)
