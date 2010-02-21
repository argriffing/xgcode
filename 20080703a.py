"""Given a newick tree, make a perturbed distance matrix.

The perturbed distance matrix will be symmetric and non-negative.
To generate the perturbed distance matrix,
distance_i will be multiplied by exp(X_i)
where X_i is a normally distributed random variable
with mean zero and standard deviation equal to the perturbation strength.
"""

from StringIO import StringIO
import random
import math

from SnippetUtil import HandlingError
import Util
import MatrixUtil
import NewickIO
import FelTree
from Form import CheckItem
import Form

def get_form():
    """
    @return: the body of a form
    """
    # define the default tree string and ordered tip labels
    tree_string = '(a:1, (b:2, d:5):1, c:4);'
    tree = NewickIO.parse(tree_string, FelTree.NewickTree)
    formatted_tree_string = NewickIO.get_narrow_newick_string(tree, 60)
    labels = list(sorted(tip.name for tip in tree.gen_tips()))
    # define the form objects
    form_objects = [
            Form.MultiLine('tree', 'newick tree',
                formatted_tree_string),
            Form.MultiLine('inlabels', 'ordered labels',
                '\n'.join(labels)),
            Form.Float('strength', 'perturbation strength',
                0.1, low_inclusive=0),
            Form.CheckGroup('options', 'output options', [
                CheckItem('perturbed', 'a perturbed distance matrix', True),
                CheckItem('distance', 'the original distance matrix'),
                CheckItem('outlabels', 'ordered labels')])]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # get the tree
    tree = NewickIO.parse(fs.tree, FelTree.NewickTree)
    alphabetically_ordered_states = list(sorted(node.name
        for node in tree.gen_tips()))
    n = len(alphabetically_ordered_states)
    if n < 2:
        raise HandlingError('the newick tree should have at least two leaves')
    # read the ordered labels
    states = Util.get_stripped_lines(StringIO(fs.inlabels))
    if len(states) > 1:
        if set(states) != set(alphabetically_ordered_states):
            msg_a = 'if ordered labels are provided, '
            msg_b = 'each should correspond to a leaf of the newick tree'
            raise HandlingError(msg_a + msg_b)
    else:
        states = alphabetically_ordered_states
    # start to prepare the reponse
    out = StringIO()
    # create the distance matrix
    D = tree.get_distance_matrix(states)
    # create the perturbed distance matrix if necessary
    if fs.strength:
        P = [row[:] for row in D]
        for i in range(n):
            for j in range(i):
                x = random.normalvariate(0, fs.strength)
                new_distance = D[i][j] * math.exp(x)
                P[i][j] = new_distance
                P[j][i] = new_distance
    else:
        P = D
    # start collecting the paragraphs
    paragraphs = []
    # show the distance matrix if requested
    if fs.perturbed:
        paragraph = StringIO()
        print >> paragraph, 'a perturbed distance matrix:'
        print >> paragraph, MatrixUtil.m_to_string(P)
        paragraphs.append(paragraph.getvalue().strip())
    # show the distance matrix if requested
    if fs.distance:
        paragraph = StringIO()
        print >> paragraph, 'the original distance matrix:'
        print >> paragraph, MatrixUtil.m_to_string(D)
        paragraphs.append(paragraph.getvalue().strip())
    # show the ordered labels if requested
    if fs.outlabels:
        paragraph = StringIO()
        print >> paragraph, 'ordered labels:'
        print >> paragraph, '\n'.join(states)
        paragraphs.append(paragraph.getvalue().strip())
    print >> out, '\n\n'.join(paragraphs)
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
