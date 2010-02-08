"""Calculate properties of the branch defined by a split of a distance matrix.

The distance matrix should represent geodesic pairwise distances
between tips of a bifurcating tree with positive branch lengths.
"""

from StringIO import StringIO

import numpy as np

from SnippetUtil import HandlingError
import SnippetUtil
import Util
import MatrixUtil
import NeighborhoodJoining
import Form

def get_form():
    """
    @return: the body of a form
    """
    # define the default distance matrix
    D = np.array([
            [0, 4, 5, 7],
            [4, 0, 7, 7],
            [5, 7, 0, 10],
            [7, 7, 10, 0]])
    states = list('abcd')
    # define the ordered and selected labels
    ordered_label_string = '\n'.join(states)
    selected_label_string = '\n'.join(['a', 'c'])
    # define the sequence of form objects
    form_objects = [
            Form.Matrix('matrix', 'distance matrix',
                D, MatrixUtil.assert_predistance),
            Form.MultiLine('labels', 'ordered labels',
                ordered_label_string),
            Form.MultiLine('selection', 'selected labels',
                selected_label_string)]
    # return the sequence of form objects
    return form_objects

def get_split_branch_arbitrary_distances(D, selection, complement):
    """
    Get distances to a translation of a point on the branch defined by the split.
    @param D: a distance matrix
    @param selection: one of the clusters, represented by a set of indices
    @param complement: the other cluser
    @return: the distance from each vertex to a point related to the branch defined by the split
    """
    # get the number of vertices
    n = len(selection | complement)
    # get the expected distance to a point in the other cluster
    E = [0]*n
    for A, B in ((selection, complement), (complement, selection)):
        for i in A:
            E[i] = sum(D[i][j] for j in B)/len(B)
    # get the mean of E
    E_bar = sum(E)/n
    # get the vector of distances to a virtual point that may or may not be on the tree
    v = [0]*n
    for A, B in ((selection, complement), (complement, selection)):
        for i in A:
            v[i] = E[i] - len(A)*E_bar/n
    return v

def get_split_branch_length(D, selection, complement):
    """
    Get the length of the branch defined by the bipartition.
    @param D: a distance matrix
    @param selection: one of the clusters, represented by a set of indices
    @param complement: the other cluser
    @return: the distance from each vertex to the root of the other subtree
    """
    # get the vector of distances to a virtual point that may or may not be on the tree
    v = get_split_branch_arbitrary_distances(D, selection, complement)
    # get two quantities that should sum to twice the length of the branch defined by the split
    a = min(v[i] + v[j] - D[i][j] for i in selection for j in selection)
    b = min(v[i] + v[j] - D[i][j] for i in complement for j in complement)
    branch_length = (a + b) / 2
    # return the branch length
    return branch_length

def get_split_branch_midpoint_distances(D, selection, complement):
    """
    Get distances to the midpoint of the the branch defined by the split.
    @param D: a distance matrix
    @param selection: one of the clusters, represented by a set of indices
    @param complement: the other cluser
    @return: the distance from each vertex to the midpoint of the branch defined by the split
    """
    n = len(selection | complement)
    # get the vector of distances to a virtual point that may or may not be on the tree
    v = get_split_branch_arbitrary_distances(D, selection, complement)
    # get two quantities that should sum to twice the length of the branch defined by the split
    a = min(v[i] + v[j] - D[i][j] for i in selection for j in selection)
    b = min(v[i] + v[j] - D[i][j] for i in complement for j in complement)
    # get the distance to the midpoint of the split
    nv = [0]*n
    for i in selection:
        nv[i] = v[i] + (b - a)/4
    for i in complement:
        nv[i] = v[i] + (a - b)/4
    return nv

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
        raise HandlingError('no ordered labels were provided')
    if len(ordered_labels) != len(set(ordered_labels)):
        raise HandlingError('the ordered labels should be unique')
    # read the set of selected labels
    selected_labels = set(Util.get_stripped_lines(StringIO(fs.selection)))
    if not selected_labels:
        raise HandlingError('no selected labels were provided')
    weird_labels = selected_labels - set(ordered_labels)
    if weird_labels:
        raise HandlingError('some selected labels are invalid: ' + str(weird_labels))
    # assert that the size of the distance matrix is compatible with the number of ordered labels
    if len(D) != len(ordered_labels):
        raise HandlingError('one ordered label should be provided for each row of the distance matrix')
    # get the set of selected indices and its complement
    n = len(D)
    selection = set(i for i, label in enumerate(ordered_labels) if label in selected_labels)
    complement = set(range(n)) - selection
    # verify that a minimum number of nodes is in the selection and the complement
    min_vertex_count = 1
    err = HandlingError('the number of vertices in each of the selected and unselected sets should be at least ' + str(min_vertex_count))
    for current_set in (selection, complement):
        if len(current_set) < min_vertex_count:
            raise err
    # begin the response
    out = StringIO()
    # get the length of the branch length defined by the split
    branch_length = get_split_branch_length(D, selection, complement)
    print >> out, 'the length of the branch defined by the bipartition:'
    print >> out, branch_length
    print >> out
    # get the distances to the midpoint
    v = get_split_branch_midpoint_distances(D, selection, complement)
    print >> out, 'distances to the midpoint of the branch defined by the bipartition:'
    for i, d in enumerate(v):
        print >> out, ordered_labels[i], ':', d
    print >> out
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
