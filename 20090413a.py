"""Given a Laplacian matrix, Schur complement out a set of taxa.
"""

from StringIO import StringIO

import numpy as np

from SnippetUtil import HandlingError
import SnippetUtil
import Util
import MatrixUtil
import Euclid
import SchurAlgebra
import Form
import FormOut

def get_form():
    """
    @return: the body of a form
    """
    # define the default matrix
    upper_matrix = np.array([
            [0, 0, 0, 0, 0, 0, 0, 0, 1/.05, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 1/.05, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 1/.2, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1/.05, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1/.05, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1/.2],
            [0, 0, 0, 0, 0, 0, 0, 0.65245955, 0, 0.72801183, 0, -0.48410939],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0.72801183, 0, -0.48410939],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 1/0.15, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.78713969],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1/0.15],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
    labels = ['a', 'b', 'c', 'm', 'n', 'p', 'x', 'y', 'ab', 'abc', 'mn', 'mnp']
    n = len(labels)
    L = [[0]*n for i in range(n)]
    for i in range(n):
        for j in range(n):
            L[i][j] = -(upper_matrix[i][j] + upper_matrix[j][i])
    row_sums = [sum(row) for row in L]
    for i in range(n):
        L[i][i] = -row_sums[i]
    # define some default strings
    ordered_label_string = '\n'.join(labels)
    selected_label_string = '\n'.join(['ab', 'abc', 'mn', 'mnp'])
    # define the sequence of form objects
    form_objects = [
            Form.Matrix('laplacian', 'laplacian matrix',
                L, MatrixUtil.assert_symmetric),
            Form.MultiLine('labels', 'ordered taxa',
                ordered_label_string),
            Form.MultiLine('selection', 'taxa to be schur complemented out',
                selected_label_string)]
    # return the sequence of form objects
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # read the matrix
    L = fs.laplacian
    # read the ordered labels
    ordered_labels = Util.get_stripped_lines(StringIO(fs.labels))
    if not ordered_labels:
        raise HandlingError('no ordered taxa were provided')
    if len(ordered_labels) != len(set(ordered_labels)):
        raise HandlingError('the ordered taxa should be unique')
    # get the label selection and its complement
    min_selected_labels = 2
    min_unselected_labels = 1
    selected_labels = set(Util.get_stripped_lines(StringIO(fs.selection)))
    if len(selected_labels) < min_selected_labels:
        raise HandlingError('at least %d taxa should be selected to be grouped' % min_selected_labels)
    # get the set of labels in the complement
    unselected_labels = set(ordered_labels) - selected_labels
    if len(unselected_labels) < min_unselected_labels:
        raise HandlingError('at least %d taxa should remain outside the selected group' % min_unselected_labels)
    # assert that no bizarre labels were selected
    weird_labels = selected_labels - set(ordered_labels)
    if weird_labels:
        raise HandlingError('some selected taxa are invalid: ' + str(weird_labels))
    # assert that the size of the distance matrix is compatible with the number of ordered labels
    if len(L) != len(ordered_labels):
        raise HandlingError('the number of listed taxa does not match the number of rows in the distance matrix')
    # get the set of selected indices and its complement
    n = len(L)
    index_selection = set(i for i, label in enumerate(ordered_labels) if label in selected_labels)
    index_complement = set(range(n)) - index_selection
    # begin the response
    out = StringIO()
    # calculate the new laplacian matrix
    L_small = SchurAlgebra.mschur(L, index_selection)
    D_small = Euclid.laplacian_to_edm(L_small)
    # print the matrices and the labels of its rows
    print >> out, 'new laplacian matrix:'
    print >> out, MatrixUtil.m_to_string(L_small)
    print >> out
    print >> out, 'new distance matrix:'
    print >> out, MatrixUtil.m_to_string(D_small)
    print >> out
    print >> out, 'new taxon labels:'
    for index in sorted(index_complement):
        print >> out, ordered_labels[index]
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
