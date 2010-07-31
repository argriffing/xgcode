"""Given a perturbed distance matrix, see how taxon grouping affects distances.

The original distance matrix represents pairwise distances on a tree,
or more generally resistance distances on a graph.
The selected taxa are replaced by a node,
and the resistance distances to this node are calculated.
If the original resistance distances are exact and imply that the
selected taxa define a subtree then pairwise distances
among the remaining taxa should not be affected.
When a perturbed distance matrix is used,
all pairwise distances should be affected.
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
import const

g_data = const.read('20100730p')

def get_form():
    """
    @return: the body of a form
    """
    # define the default distance matrix
    # this is from figure two of a paper called why neighbor joining works
    D = MatrixUtil.read_matrix(Util.get_stripped_lines(g_data.splitlines()))
    states = list('xyabcmnp')
    # define some default strings
    ordered_label_string = '\n'.join(states)
    selected_label_string = '\n'.join(['m', 'n'])
    # define the sequence of form objects
    form_objects = [
            Form.Matrix('matrix', 'distance matrix',
                D, MatrixUtil.assert_predistance),
            Form.MultiLine('labels', 'ordered taxa',
                ordered_label_string),
            Form.MultiLine('selection', 'taxa to be grouped',
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
    D = fs.matrix
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
    if len(D) != len(ordered_labels):
        raise HandlingError('the number of listed taxa does not match the number of rows in the distance matrix')
    # get the set of selected indices and its complement
    n = len(D)
    index_selection = set(i for i, label in enumerate(ordered_labels) if label in selected_labels)
    index_complement = set(range(n)) - index_selection
    # begin the response
    out = StringIO()
    # get the ordered list of sets of indices to merge
    merged_indices = SchurAlgebra.vmerge([set([x]) for x in range(n)], index_selection)
    # calculate the new distance matrix
    L = Euclid.edm_to_laplacian(D)
    L_merged = SchurAlgebra.mmerge(L, index_selection)
    D_merged = Euclid.laplacian_to_edm(L_merged)
    # print the output distance matrix and the labels of its rows
    print >> out, 'new distance matrix:'
    print >> out, MatrixUtil.m_to_string(D_merged)
    print >> out
    print >> out, 'new taxon labels:'
    for merged_index_set in merged_indices:
        if len(merged_index_set) == 1:
            print >> out, ordered_labels[merged_index_set.pop()]
        else:
            print >> out, '{' + ', '.join(selected_labels) + '}'
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
