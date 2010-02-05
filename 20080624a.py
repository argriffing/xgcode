"""Calculate ancestral state distributions given leaf states.
"""

from StringIO import StringIO

import numpy

from SnippetUtil import HandlingError
import SnippetUtil
import Util
import MatrixUtil
import Newick
import RateMatrix
import Form

def get_form():
    """
    @return: the body of a form
    """
    # define the default newick tree
    tree_string = '((((a:.05, b:.05)L1:.15, c:.2)L2:.8, x:1)L3:.5, (((m:.05, n:.05)R1:.15, p:.2)R2:.8, y:1)R3:.5)root;'
    tree = Newick.parse(tree_string, Newick.NewickTree)
    formatted_tree_string = Newick.get_narrow_newick_string(tree, 60)
    # define the default rate matrix and its ordered states
    R = numpy.array([
        [-.3, .1, .1, .1],
        [.1, -.3, .1, .1],
        [.1, .1, -.3, .1],
        [.1, .1, .1, -.3]])
    ordered_states = list('ACGT')
    # define the default leaf state assigments
    leaf_assignment_lines = [
            'a : A',
            'b : A',
            'c : A',
            'x : A',
            'm : C', 
            'n : C', 
            'p : C', 
            'y : C']
    # define the form objects
    form_objects = [
            Form.MultiLine('tree', 'newick tree with branch lengths', formatted_tree_string),
            Form.Matrix('rate_matrix', 'rate matrix', R, MatrixUtil.assert_rate_matrix),
            Form.MultiLine('states', 'ordered states', '\n'.join(ordered_states)),
            Form.MultiLine('assignments', 'leaf states', '\n'.join(leaf_assignment_lines))]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # get the tree
    tree = Newick.parse(fs.tree, Newick.NewickTree)
    tree.assert_valid()
    # read the ordered states
    ordered_states = list(Util.stripped_lines(StringIO(fs.states)))
    # read the matrix from the form data
    R = fs.rate_matrix
    if len(R) < 2:
        raise HandlingError('the rate matrix should have at least two rows')
    if len(ordered_states) != len(R):
        raise HandlingError('the number of ordered states should be the same as the number of rows in the matrix')
    # get the dictionary mapping taxa to states
    taxon_to_state = SnippetUtil.get_generic_dictionary(StringIO(fs.assignments), 'taxon name', 'state name', ordered_states)
    # set the states for each of the tree tips
    for node in tree.gen_tips():
        node.state = taxon_to_state[node.name]
    # create the rate matrix object
    rate_matrix_object = RateMatrix.RateMatrix(R.tolist(), ordered_states)
    # repeatedly reroot and calculate root state distributions
    internal_nodes = list(tree.gen_internal_nodes())
    for node in internal_nodes:
        tree.reroot(node)
        rate_matrix_object.add_probabilities(tree)
        weights = [node.state_to_subtree_prob[state] for state in ordered_states]
        node.state_distribution = Util.weights_to_distribution(weights)
    # define the response
    out = StringIO()
    # show the ancestral state distributions
    for node in tree.gen_internal_nodes():
        if node.name:
            print >> out, node.name, ':', '\t'.join(str(p) for p in node.state_distribution)
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
