"""Sample a nucleotide alignment given a tree and a mixture model.

The mixture is scaled so that the branch lengths in the newick tree
are the expected number of substitutions on the branch.
The rows and columns of the rate matrices
are ordered alphabetically by nucleotide.
"""

from StringIO import StringIO

import numpy as np

from SnippetUtil import HandlingError
import SnippetUtil
import Newick
import PhyLikelihood
import RateMatrix
import MatrixUtil
import SubModel
import Form

g_nt_matrix_a = np.array([
    [-3, 1, 1, 1], [1, -3, 1, 1], [1, 1, -3, 1], [1, 1, 1, -3]])
g_nt_matrix_b = np.array([
    [-3, 1, 1, 1], [1, -3, 1, 1], [1, 1, -3, 1], [1, 1, 1, -3]])
g_nt_matrix_c = np.array([
    [-3, 1, 1, 1], [10, -12, 1, 1], [10, 1, -12, 1], [10, 1, 1, -12]])
g_weight_a = 1
g_weight_b = 2
g_weight_c = 3

def get_form():
    """
    @return: the body of a form
    """
    # define the tree string
    tree_string = '(((Human:0.1, Chimpanzee:0.2):0.8, Gorilla:0.3):0.7, Orangutan:0.4, Gibbon:0.5);'
    tree = Newick.parse(tree_string, Newick.NewickTree)
    formatted_tree_string = Newick.get_narrow_newick_string(tree, 60)
    # define the form objects
    form_objects = [
            Form.MultiLine('tree', 'newick tree', formatted_tree_string),
            Form.Integer('ncols', 'sample this many columns',
                100, low=1, high=1000),
            Form.Matrix('matrix_a', 'first nucleotide rate matrix',
                g_nt_matrix_a),
            Form.Float('weight_a', 'first mixture weight',
                g_weight_a, low_inclusive=0),
            Form.Matrix('matrix_b', 'second nucleotide rate matrix',
                g_nt_matrix_b),
            Form.Float('weight_b', 'second mixture weight',
                g_weight_b, low_inclusive=0),
            Form.Matrix('matrix_c', 'third nucleotide rate matrix',
                g_nt_matrix_c),
            Form.Float('weight_c', 'third mixture weight',
                g_weight_c, low_inclusive=0)]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # get the tree
    tree = Newick.parse(fs.tree, Newick.NewickTree)
    tree.assert_valid()
    # get the mixture weights
    weights = [fs.weight_a, fs.weight_b, fs.weight_c]
    # get the matrices
    matrices = [fs.matrix_a, fs.matrix_b, fs.matrix_c]
    for R in matrices:
        if R.shape != (4,4):
            msg = 'expected each nucleotide rate matrix to be 4x4'
            raise HandlingError(msg)
    # create the mixture proportions
    weight_sum = sum(weights)
    mixture_proportions = [weight / weight_sum for weight in weights]
    # create the rate matrix objects
    ordered_states = list('ACGT')
    rate_matrix_objects = []
    for R in matrices:
        rate_matrix_object = RateMatrix.RateMatrix(R.tolist(), ordered_states)
        rate_matrix_objects.append(rate_matrix_object)
    # create the mixture model
    mixture_model = SubModel.MixtureModel(mixture_proportions,
            rate_matrix_objects)
    # normalize the mixture model
    mixture_model.normalize()
    # simulate the alignment
    try:
        alignment = PhyLikelihood.simulate_alignment(tree,
                mixture_model, fs.ncols)
    except PhyLikelihood.SimulationError, e:
        raise HandlingError(e)
    # get the alignment string
    arr = []
    for node in tree.gen_tips():
        arr.append(alignment.get_fasta_sequence(node.name))
    alignment_string = '\n'.join(arr)
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, alignment_string
