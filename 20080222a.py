"""Visualize a matrix as a heat map.
"""

from StringIO import StringIO

import numpy

from SnippetUtil import HandlingError
import SnippetUtil
import EnglishModel
import HeatMap
import MatrixUtil
import Form

def get_form():
    """
    @return: the body of a form
    """
    # define the default matrix lines
    dictionary_rate_matrix = EnglishModel.get_transition_matrix()
    labels = list(sorted(set(a for a, b in dictionary_rate_matrix)))
    T = numpy.array(MatrixUtil.dict_to_row_major(dictionary_rate_matrix, labels, labels))
    # define the form objects
    form_objects = [
            Form.Matrix('matrix', 'matrix', T),
            Form.Integer('maxcategories', 'maximum number of categories', 5, low=2)]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # create the html response representing the heat map
    M = fs.matrix
    heatmap = HeatMap.HeatMap(M.tolist(), fs.maxcategories)
    renderer = HeatMap.PreHeatMap(heatmap)
    html_string = renderer.get_example_html()
    # return the response
    response_headers = [('Content-Type', 'text/html')]
    return response_headers, html_string
