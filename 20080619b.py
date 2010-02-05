"""Calculate an objective function of a given graph cut.
"""

from StringIO import StringIO

from scipy import linalg
import numpy

from SnippetUtil import HandlingError
import Util
import MatrixUtil
import Form

def get_form():
    """
    @return: the body of a form
    """
    # define the default matrix, the ordered labels, and the selected labels
    M = numpy.array([
        [0, 4, 5, 7],
        [4, 0, 7, 7],
        [5, 7, 0, 10],
        [7, 7, 10, 0]])
    labels = list('abcd')
    selection = list('ac')
    # define the form objects
    form_objects = [
            Form.Matrix('matrix', 'weight matrix', M, MatrixUtil.assert_symmetric),
            Form.MultiLine('labels', 'ordered_labels', '\n'.join(labels)),
            Form.MultiLine('selection', 'selected labels', '\n'.join(selection)),
            Form.RadioGroup('objective', 'bipartition objective function', [
                Form.RadioItem('min', 'min cut', True),
                Form.RadioItem('conductance', 'min conductance cut')])]
    return form_objects

def get_conductance(selection, affinity):
    """
    @param selection: a cluster defined as a set of indices
    @param affinity: a row major affinity matrix
    """
    # count the states
    n = len(affinity)
    # define the assignment vector
    assignment = []
    for i in range(n):
        if i in selection:
            assignment.append(1)
        else:
            assignment.append(-1)
    # get the affinity weights
    between_clusters = 0.0
    adjoining_cluster_a = 0.0
    adjoining_cluster_b = 0.0
    for i in range(n):
        for j in range(n):
            if i < j:
                weight = affinity[i][j]
                if assignment[i] != assignment[j]:
                    between_clusters += weight
                    adjoining_cluster_a += weight
                    adjoining_cluster_b += weight
                elif assignment[i] == 1:
                    adjoining_cluster_a += weight
                elif assignment[i] == -1:
                    adjoining_cluster_b += weight
    return between_clusters / min(adjoining_cluster_a, adjoining_cluster_b)

def get_cut(selection, affinity):
    """
    @param selection: a cluster defined as a set of indices
    @param affinity: a row major affinity matrix
    """
    n = len(affinity)
    complement = set(range(n)) - selection
    return sum(affinity[i][j] for i in selection for j in complement)

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # read the matrix
    M = fs.matrix
    # read the ordered labels
    ordered_labels = list(Util.stripped_lines(StringIO(fs.labels)))
    # read the set of selected labels
    selected_labels = set(Util.stripped_lines(StringIO(fs.selection)))
    # get the set of selected indices
    selection = set(i for i, label in enumerate(ordered_labels) if label in selected_labels)
    # get the value of the objective function
    if fs.min:
        value = get_cut(selection, M.tolist())
    elif fs.conductance:
        value = get_conductance(selection, M.tolist())
    # start to prepare the reponse
    out = StringIO()
    # show the value of the objective function
    print >> out, 'objective function value:'
    print >> out, value
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
