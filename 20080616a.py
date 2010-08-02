"""Given an adjacency matrix, find an optimal bipartition.

Given a weighted adjacency matrix,
find the optimal bipartition using some objective function.
"""

from StringIO import StringIO

import numpy as np

from SnippetUtil import HandlingError
import Util
import MatrixUtil
import StoerWagner
import Clustering
import Form
import FormOut

def get_form():
    """
    @return: the body of a form
    """
    # define the default matrix and its ordered labels
    A = np.array(StoerWagner.g_stoer_wagner_affinity)
    ordered_labels = [str(i+1) for i in range(len(A))]
    # define the form objects
    form_objects = [
            Form.Matrix('matrix', 'weighted adjacency matrix',
                A, MatrixUtil.assert_weighted_adjacency),
            Form.MultiLine('labels', 'ordered labels',
                '\n'.join(ordered_labels)),
            Form.RadioGroup('objective', 'bipartition objective function', [
                Form.RadioItem('min', 'min cut', True),
                Form.RadioItem('conductance', 'min conductance cut')])]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):
    # read the weighted adjacency matrix
    A = fs.matrix
    # read the labels
    ordered_labels = Util.get_stripped_lines(StringIO(fs.labels))
    # Assert that the number of labels
    # is compatible with the shape of the matrix.
    n = len(A)
    if len(ordered_labels) != n:
        msg = 'the number of labels does not match the size of the matrix'
        raise HandlingError(msg)
    # get the best objective function value and the corresponding best cluster
    if fs.conductance:
        max_size = 20
        if n > max_size:
            msg_a = 'for the min conductance objective function please '
            msg_b = 'limit the size of the matrix to %d rows' % max_size
            raise HandlingError(msg_a + msg_b)
        pairs = [(get_conductance(assignment, A), assignment)
                for assignment in Clustering.gen_assignments(n)]
        best_objective, best_assignment = min(pairs)
        best_cluster = set(i for i in range(n) if best_assignment[i] == 1)
    if fs.min:
        best_cluster = StoerWagner.stoer_wagner_min_cut(A.tolist())
        complement = set(range(n)) - best_cluster
        best_objective = sum(A[i][j] for i in best_cluster for j in complement)
    # get the smaller of the two clusters
    complement = set(range(n)) - best_cluster
    small_cluster = min((len(best_cluster), best_cluster),
            (len(complement), complement))[1]
    # start to prepare the reponse
    out = StringIO()
    print >> out, 'smallest cluster defined by the bipartition:'
    for index in sorted(small_cluster):
        print >> out, ordered_labels[index]
    print >> out, ''
    print >> out, 'objective function value:'
    print >> out, best_objective
    # write the response
    return out.getvalue()

def get_conductance(assignment, affinity):
    """
    The assignment vector has elements equal to 1 or -1.
    These elements define the cluster membership.
    @param assignment: defines the cluster membership
    @param affinity: an affinity matrix
    """
    # if the assignment defines a trivial partition then return None
    if set(assignment) != set([-1, 1]):
        return None
    # count the states
    n = len(assignment)
    # get the affinity counts
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
