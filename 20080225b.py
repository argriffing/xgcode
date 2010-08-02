"""Given a rate matrix, get the stationary distribution.
"""

import numpy as np

from SnippetUtil import HandlingError
import MatrixUtil
import RateMatrix
import Form
import FormOut

#FIXME numpy may not be necessary

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

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):
    # read the matrix from the form data
    R = fs.matrix
    # get the stationary distribution of the rate matrix
    try:
        v = RateMatrix.get_stationary_distribution(R.tolist())
    except RateMatrix.RateMatrixError as e:
        msg = 'error calculating the stationary distribution: ' + str(e)
        raise HandlingError(msg)
    # return the stationary distribution string
    return '\n'.join(str(x) for x in v) + '\n'
