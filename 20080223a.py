"""Visualize a codon rate matrix as a heat map.
"""

from StringIO import StringIO

import numpy as np

from SnippetUtil import HandlingError
import SnippetUtil
import MatrixUtil
import RateMatrix
import HeatMap
import Codon
import Form

def get_form():
    """
    @return: the body of a form
    """
    # define the default rate matrix
    dictionary_rate_matrix = RateMatrix.get_sample_codon_rate_matrix()
    labels = Codon.g_sorted_non_stop_codons
    R = MatrixUtil.dict_to_row_major(dictionary_rate_matrix, labels, labels)
    R = np.array(R)
    # define the form objects
    form_objects = [
            Form.Matrix('matrix', 'codon rate matrix',
                R, MatrixUtil.assert_rate_matrix),
            Form.Integer('maxcategories', 'maximum number of categories',
                5, low=2)]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # read the matrix from the form data
    R = fs.matrix
    nrows, ncols = R.shape
    # assert that the number of rows and columns is valid for a codon matrix
    states = Codon.g_sorted_non_stop_codons
    if nrows != len(states):
        msg = 'expected %d rows but got %d' % (len(states), nrows)
        raise HandlingError(msg)
    if ncols != len(states):
        msg = 'expected %d columns but got %d' % (len(states), ncols)
        raise HandlingError(msg)
    # define the row and column labels
    labels = []
    for codon in states:
        label = '%s.%s.' % (Codon.g_codon_to_aa_letter[codon], codon)
        labels.append(label)
    row_labels = labels
    column_labels = labels
    # initialize the base class with this row major matrix
    heatmap = HeatMap.LabeledHeatMap(R.tolist(), fs.maxcategories,
            row_labels, column_labels)
    renderer = HeatMap.PreHeatMap(heatmap)
    html_string = renderer.get_example_html()
    # return the response
    response_headers = [('Content-Type', 'text/html')]
    return response_headers, html_string
