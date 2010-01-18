"""Analyze properties of the augmented distance matrix of a tree.

The augmented distance matrix includes internal nodes.
The M matrix is a function of the inverted distance matrix.
NOTE: Here the term 'M matrix' is not necessarily used in its technical sense.
I should probably change the notation.
"""

import StringIO

import numpy
import scipy.linalg as linalg

from SnippetUtil import HandlingError
import MatrixUtil
import Clustering
import NewickIO
import FelTree
import Form

g_sample_tree_string = """(
    ((a:0.05, b:0.05)AB:0.15, c:0.2)ABC:0.8,
    x:1,
    (((m:0.05, n:0.05)MN:0.15, p:0.2)MNP:0.8, y:1)MNPY:1
)ABCX;
"""

def get_form():
    """
    @return: the body of a form
    """
    # define the default tree string
    tree = NewickIO.parse(g_sample_tree_string, FelTree.NewickTree)
    formatted_tree_string = NewickIO.get_narrow_newick_string(tree, 60)
    # define the form objects
    form_objects = [
            Form.MultiLine('tree', 'newick tree', formatted_tree_string),
            Form.RadioGroup('matrix', 'use these nodes for the distance matrix', [
                Form.RadioItem('standard', 'tips only', True),
                Form.RadioItem('augmented', 'all nodes'),
                Form.RadioItem('named', 'all named nodes')]),
            Form.CheckGroup('output_options' , 'output options', [
                Form.CheckItem('show_split', 'exact criterion partition', True),
                Form.CheckItem('show_value', 'exact criterion value', True),
                Form.CheckItem('show_value_minus_trace', 'exact criterion value minus trace', True),
                Form.CheckItem('show_fiedler_split', 'show the spectral sign partition', True),
                Form.CheckItem('show_fiedler_eigenvector', 'show the eigenvector of interest', True),
                Form.CheckItem('show_labels', 'ordered labels', True),
                Form.CheckItem('show_distance_matrix', 'distance matrix', True),
                Form.CheckItem('show_M_matrix', 'M matrix', True)])]
    return form_objects

def get_eigenvector_of_interest(row_major_matrix):
    """
    This gets a left eigenvector because of the standard format of the matrix.
    @param row_major_matrix: a matrix
    @return: an eigenvector
    """
    R = numpy.array(row_major_matrix)
    w, vl, vr = linalg.eig(R, left=True, right=True)
    eigenvalue_info = list(sorted((abs(x), i) for i, x in enumerate(w)))
    stationary_eigenvector_index = eigenvalue_info[0][1]
    first_axis_eigenvector_index = eigenvalue_info[1][1]
    second_axis_eigenvector_index = eigenvalue_info[2][1]
    return vl.T[first_axis_eigenvector_index]

def set_to_string(set_in):
    """
    Use a less verbose way to show a set as a string.
    @param set_in: a set
    @return: a string that shows the set
    """
    return '{' + ', '.join(set_in) + '}'

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # start writing the response type
    response_headers = []
    # read the tree
    tree = NewickIO.parse(fs.tree, FelTree.NewickTree)
    # create a putative list of nodes
    putative_nodes = []
    putative_nodes.extend(list(tree.gen_tips()))
    if not fs.standard:
        putative_nodes.extend(list(tree.gen_internal_nodes()))
    # get the ordered ids and ordered names of the selected nodes in the tree
    ordered_name_id_pairs = []
    for node in putative_nodes:
        name = node.get_name()
        if fs.named and not name:
            continue
        ordered_name_id_pairs.append((name, id(node)))
    ordered_ids = [id_ for name, id_ in ordered_name_id_pairs]
    ordered_names = [name for name, id_ in ordered_name_id_pairs]
    id_to_index = dict((id_, i) for i, id_ in enumerate(ordered_ids))
    # assert that names are non-empty
    for name in ordered_names:
        if not name:
            raise HandlingError('each node must be named')
    # assert that names are unique
    n = len(ordered_ids)
    if len(set(ordered_names)) != n:
        raise HandlingError('each node must be uniquely named')
    # get the R matrix from the tree; this is -1/2 times the laplacian matrix
    if fs.standard:
        D = tree.get_distance_matrix(ordered_names)
    elif fs.augmented:
        D = tree.get_full_distance_matrix(ordered_ids)
    elif fs.named:
        D = tree.get_partial_distance_matrix(ordered_ids)
    R = Clustering.get_R_balaji(D)
    R_trace = sum(R[i][i] for i in range(n))
    # get the best partition and partition value
    value_Y_pairs = []
    for Y in Clustering.gen_assignments(n):
        value = Clustering.get_exact_criterion(R, Y)
        value_Y_pairs.append((value, Y))
    best_value, best_Y = max(value_Y_pairs)
    # convert the best Y vector to a partition
    pos_set = set(ordered_names[i] for i, elem in enumerate(best_Y) if elem > 0)
    neg_set = set(ordered_names[i] for i, elem in enumerate(best_Y) if elem < 0)
    # get fiedler split information
    fiedler_eigenvector = get_eigenvector_of_interest(R)
    fiedler_pos_set = set(ordered_names[i] for i, elem in enumerate(fiedler_eigenvector) if elem > 0)
    fiedler_neg_set = set(ordered_names[i] for i, elem in enumerate(fiedler_eigenvector) if elem < 0)
    # write the paragraphs
    paragraphs = []
    if fs.show_split:
        lines = [
                'exact criterion partition:',
                str(list(best_Y)),
                set_to_string((set_to_string(neg_set), set_to_string(pos_set)))
                ]
        paragraphs.append('\n'.join(lines))
    if fs.show_value:
        lines = [
                'exact criterion value:',
                str(best_value)
                ]
        paragraphs.append('\n'.join(lines))
    if fs.show_value_minus_trace:
        lines = [
                'exact criterion value minus trace:',
                str(best_value - R_trace)
                ]
        paragraphs.append('\n'.join(lines))
    if fs.show_fiedler_split:
        lines = [
                'spectral sign partition:',
                set_to_string((set_to_string(fiedler_neg_set), set_to_string(fiedler_pos_set)))
                ]
        paragraphs.append('\n'.join(lines))
    if fs.show_fiedler_eigenvector:
        lines = [
                'eigenvector of interest:',
                str(list(fiedler_eigenvector))
                ]
        paragraphs.append('\n'.join(lines))
    if fs.show_labels:
        lines = ['ordered labels:'] + ordered_names
        paragraphs.append('\n'.join(lines))
    if fs.show_distance_matrix:
        if fs.augmented:
            title = 'augmented distance matrix:'
        elif fs.standard:
            title = 'distance matrix:'
        elif fs.named:
            title = 'distance matrix:'
        lines = [
                title,
                MatrixUtil.m_to_string(D)
                ]
        paragraphs.append('\n'.join(lines))
    if fs.show_M_matrix:
        lines = [
                'M matrix:',
                MatrixUtil.m_to_string(R)
                ]
        paragraphs.append('\n'.join(lines))
    # write the response
    out = StringIO.StringIO()
    print >> out, '\n\n'.join(paragraphs)
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
