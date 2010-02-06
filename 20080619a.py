"""Given a distance matrix, get the branch length between two subtrees.

Given a distance matrix, find the length of the branch connecting two subtrees.
This procedure is based on R code by Eric Stone.
"""

from StringIO import StringIO

import numpy as np

from SnippetUtil import HandlingError
import Util
import MatrixUtil
import Form

def get_form():
    """
    @return: the body of a form
    """
    # Define the default distance matrix, the ordered labels,
    # and the selected labels.
    D = np.array([
        [0, 4, 5, 7],
        [4, 0, 7, 7],
        [5, 7, 0, 10],
        [7, 7, 10, 0]])
    ordered_labels = list('abcd')
    selected_labels = list('ac')
    # define the form objects
    form_objects = [
            Form.Matrix('matrix', 'distance matrix',
                D, MatrixUtil.assert_predistance),
            Form.MultiLine('labels', 'ordered labels',
                '\n'.join(ordered_labels)),
            Form.MultiLine('selection', 'selected labels',
                '\n'.join(selected_labels))]
    return form_objects

def get_branch_length_a(D, selection, complement):
    """
    This is how Eric Stone finds the length of a branch joining two subtrees.
    @param D: a distance matrix
    @param selection: one of the clusters, represented by a set of indices
    @param complement: the other cluser
    @return: the length of the branch that joins the clusters
    """
    # get the number of vertices
    n = len(selection | complement)
    # Get the vector of distances to a virtual point
    # which may or may not be on the tree.
    cut_distance = sum(D[i][j] for i in selection for j in complement)
    v = [0]*n
    for A, B in ((selection, complement), (complement, selection)):
        for i in A:
            partial_cut_distance = sum(D[i][j] for j in B)
            v[i] = (n*partial_cut_distance - cut_distance) / (n*len(B))
    # get the length of the branch connecting the two leaf subsets
    branch_length = 0
    for A in (selection, complement):
        branch_length += min(v[i] + v[j] - D[i][j] for i in A for j in A)/2
    return branch_length

def get_branch_length_b(D, selection, complement):
    """
    This is another way to find the length of a branch joining two subtrees.
    @param D: a distance matrix
    @param selection: one of the clusters, represented by a set of indices
    @param complement: the other cluser
    @return: the length of the branch that joins the clusters
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
    # Get the vector of distances to a virtual point
    # which may or may not be on the tree.
    v = [0]*n
    for A, B in ((selection, complement), (complement, selection)):
        for i in A:
            v[i] = E[i] - len(A)*E_bar/n
    # get the branch length
    branch_length = 0
    for A in (selection, complement):
        branch_length += min(v[i] + v[j] - D[i][j] for i in A for j in A)/2
    return branch_length

def get_branch_length_c(D, selection, complement):
    """
    This is a third way to find the length of a branch joining two subtrees.
    The idea is to define three equations and three unknowns.
    The three unknowns are:
        - the expected distance from a node in the selection subtree
          to the subtree root
        - the expected distance from a node in the complement subtree
          to the subtree root
        - the branch length connecting the selection subtree
          to the complement subtree
    @param D: a distance matrix
    @param selection: one of the clusters, represented by a set of indices
    @param complement: the other cluser
    @return: the length of the branch that joins the clusters
    """
    # get the number of vertices
    n = len(selection | complement)
    # get the expected distance to a point in the other subtree
    E = [0]*n
    for A, B in ((selection, complement), (complement, selection)):
        for i in A:
            E[i] = sum(D[i][j] for j in B) / len(B)
    # Get the expected distance from a point in one subtree
    # to a point in the other subtree.
    alpha = sum(E) / len(E)
    # Get the expected distance from a point in the selection subtree
    # to the root of the other subtree.
    beta = min(E[i] + E[j] - D[i][j]
            for i in complement for j in complement) / 2
    # Get the expected distance from a point in the complement subtree
    # to the root of the other subtree.
    gamma = min(E[i] + E[j] - D[i][j]
            for i in selection for j in selection) / 2
    # Using the three equations alpha, beta, and gamma,
    # solve for the three variables.
    branch_length = beta + gamma - alpha
    expected_selection_distance = beta - branch_length
    expected_complement_distance = gamma - branch_length
    # return the branch length
    return branch_length


def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # read the matrix
    D = fs.matrix
    # read the ordered labels
    ordered_labels = Util.get_stripped_lines(StringIO(fs.labels))
    # read the set of selected labels
    selected_labels = Util.get_stripped_lines(StringIO(fs.selection))
    # start to prepare the reponse
    out = StringIO()
    # get the set of selected indices and its complement
    n = len(D)
    selection = set(i
            for i, x in enumerate(ordered_labels) if x in selected_labels)
    complement = set(range(n)) - selection
    branch_length = get_branch_length_c(D, selection, complement)
    print >> out, 'branch length:'
    print >> out, branch_length
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
