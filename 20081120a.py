"""Compare part of a tree Fiedler vector to that of a related graph.

Compare a subvector of the Fiedler vector of a tree
to the Fiedler vector of a related graph.
"""

from StringIO import StringIO

from scipy import linalg
import numpy as np

from SnippetUtil import HandlingError
import SnippetUtil
import Form
import FormOut
import MatrixUtil
import Clustering
import NewickIO
import FelTree
import const

g_default_string = const.read('20100730q')

def get_form():
    """
    @return: a list of form objects
    """
    tree = NewickIO.parse(g_default_string, FelTree.NewickTree)
    formatted_tree_string = NewickIO.get_narrow_newick_string(tree, 60)
    return [Form.MultiLine('tree', 'tree', formatted_tree_string)]

def get_form_out():
    return FormOut.Report()

def get_eigenvectors(row_major_matrix):
    """
    Return two eigenvectors.
    This gets a couple of left eigenvectors
    because of the standard format of rate matrices.
    @param row_major_matrix: this is supposed to be a rate matrix
    @return: a pair of eigenvectors
    """
    R = np.array(row_major_matrix)
    w, vl, vr = linalg.eig(R, left=True, right=True)
    eigenvalue_info = list(sorted((abs(x), i) for i, x in enumerate(w)))
    stationary_eigenvector_index = eigenvalue_info[0][1]
    first_eigenvector_index = eigenvalue_info[1][1]
    second_eigenvector_index = eigenvalue_info[2][1]
    return vl.T[first_eigenvector_index], vl.T[second_eigenvector_index]

def get_response(fs):
    """
    @param fs: a FieldStorage object decorated with field values
    @return: a (response_headers, response_text) pair
    """
    # start writing the response type
    response_headers = []
    # read the tree
    tree = NewickIO.parse(fs.tree, FelTree.NewickTree)
    # get ordered identifiers
    ordered_tip_name_id_pairs = list(sorted(set((node.get_name(), id(node))
        for node in tree.gen_tips())))
    ordered_tip_names, ordered_tip_ids = zip(*ordered_tip_name_id_pairs)
    ordered_internal_ids = [id(node)
            for node in tree.preorder() if not node.is_tip()]
    ordered_ids = list(ordered_tip_ids) + ordered_internal_ids
    # get the distance matrices
    full_D = tree.get_partial_distance_matrix(ordered_ids)
    partial_D = tree.get_partial_distance_matrix(ordered_tip_ids)
    # get the balaji matrices
    full_R = Clustering.get_R_balaji(full_D)
    partial_R = Clustering.get_R_balaji(partial_D)
    # Get the fiedler eigenvector and another eigenvector
    # for the full and the partial balaji matrices.
    full_va, full_vb = get_eigenvectors(full_R)
    partial_va, partial_vb = get_eigenvectors(partial_R)
    # create the response
    out = StringIO()
    print >> out, 'Fiedler vector associated with the graph'
    print >> out, 'for which the internal nodes are hidden:'
    print >> out, str(tuple(partial_va))
    print >> out
    print >> out, 'The tip subvector of the Fiedler vector'
    print >> out, 'associated with the graph of the full tree:'
    print >> out, str(tuple(full_va[:len(ordered_tip_ids)]))
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
