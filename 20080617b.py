"""Given an adjacency matrix, use its eigendecomposition to find a bipartition.

Given a weighted adjacency matrix,
use its eigendecomposition to find a bipartition.
"""

from StringIO import StringIO

import numpy as np

from SnippetUtil import HandlingError
import Util
import MatrixUtil
import Form
import FormOut

def get_form():
    """
    @return: the body of a form
    """
    # define the default matrix and its ordered labels
    A = np.array([
            [0, 2, 2],
            [2, 0, 2],
            [2, 2, 0]])
    labels = list('abc')
    # define the form objects
    form_objects = [
            Form.Matrix('matrix', 'weighted adjacency matrix',
                A, MatrixUtil.assert_weighted_adjacency),
            Form.MultiLine('labels', 'ordered labels', '\n'.join(labels)),
            Form.RadioGroup('criterion', 'threshold criterion', [
                Form.RadioItem('sign', 'sign cut', True),
                Form.RadioItem('median', 'median cut')])]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):
    # read the matrix from the form data
    A = fs.matrix
    n = len(A)
    if n < 3:
        raise HandlingError('expected at least a 3x3 matrix')
    # read the ordered labels
    ordered_labels = Util.get_stripped_lines(StringIO(fs.labels))
    # do the eigendecomposition
    w, v = np.linalg.eigh(A)
    eigenvalue_info = list(sorted((abs(x), i) for i, x in enumerate(w)))
    stationary_eigenvector_index = eigenvalue_info[0][1]
    fiedler_eigenvector_index = eigenvalue_info[1][1]
    fiedler_eigenvector = v.T[fiedler_eigenvector_index]
    # respond according to the criterion
    out = StringIO()
    selected_indices = None
    if fs.sign:
        selected_indices = set(i
                for i, el in enumerate(fiedler_eigenvector) if el < 0)
    elif fs.median:
        element_index_pairs = list(sorted((el, i)
            for i, el in enumerate(fiedler_eigenvector)))
        selected_indices = set(index
                for el, index in element_index_pairs[:n/2])
    if not selected_indices:
        raise HandlingError('found a degenerate bipartition')
    # show the bipartition
    selection = set(ordered_labels[i] for i in selected_indices)
    complement = set(ordered_labels) - selection
    smallest_cluster = min((len(selection), selection),
            (len(complement), complement))[1]
    print >> out, 'labels belonging to the smaller cluster:'
    for label in sorted(smallest_cluster):
        print >> out, label
    # return the response
    return out.getvalue()
