"""Compute some splits of a tree.
"""

from StringIO import StringIO

import numpy as np

from SnippetUtil import HandlingError
import MatrixUtil
import Form
import FormOut
import NewickIO
import FelTree
import Euclid
import BuildTreeTopology
import SchurAlgebra

def get_form():
    """
    @return: the body of a form
    """
    # Define the default tree string with branch lengths
    # and named internal nodes.
    tree_string = '(a:2, (b:2, c:9)g:4, ((d:1, e:3)i:7, f:2)j:1)h;'
    tree = NewickIO.parse(tree_string, FelTree.NewickTree)
    formatted_tree_string = NewickIO.get_narrow_newick_string(tree, 60)
    # define the form objects
    form_objects = [
            Form.MultiLine('tree',
                'newick tree with branch lengths', formatted_tree_string),
            Form.Float('scaling_factor',
                'show Laplacian matrices scaled by this amount', 13320),
            Form.RadioGroup('matrix_format', 'matrix format', [
                Form.RadioItem('plain_matrix', 'plain', True),
                Form.RadioItem('latex_matrix', 'LaTeX')])]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_reciprocal_matrix(M):
    """
    @param M: some matrix
    @return: a new matrix where each element is inverted unless it is zero
    """
    arr = []
    for line in M:
        row = [0 if not value else 1/value for value in line]
        arr.append(row)
    return np.array(arr)

def get_full_tree_message(tree, m_to_string):
    """
    In this function we find the Fiedler split of the full tree.
    @param tree: each node in this tree must have a name
    @param m_to_string: a function that converts a matrix to a string
    @return: a message about the split of the tips of the tree induced by the fiedler vector
    """
    out = StringIO()
    # get the alphabetically ordered names
    ordered_names = list(sorted(node.get_name() for node in tree.preorder()))
    # get the corresponding ordered ids
    name_to_id = dict((node.get_name(), id(node)) for node in tree.preorder())
    ordered_ids = [name_to_id[name] for name in ordered_names]
    # get the full weighted adjacency matrix
    A = np.array(tree.get_affinity_matrix(ordered_ids))
    print >> out, 'the weighted reciprocal adjacency matrix of the full tree:'
    print >> out, m_to_string(get_reciprocal_matrix(A))
    print >> out
    # get the full Laplacian matrix
    L = Euclid.adjacency_to_laplacian(A)
    # get the fiedler split
    v = BuildTreeTopology.laplacian_to_fiedler(L)
    print >> out, 'the Fiedler split of the full tree:'
    for name, value in zip(ordered_names, v):
        print >> out, name, ':', value
    return out.getvalue().strip()

def label_set_to_string(label_set, label_to_name):
    """
    @param label_set: a set of one or more integers
    @param label_to_name: the names associated with the labels
    """
    if len(label_set) == 1:
        label, = label_set
        return label_to_name[label]
    else:
        names = [label_to_name[label] for label in sorted(label_set)]
        return '{' + ', '.join(names) + '}'

def get_child_messages(L, eigensplit, ordered_tip_names, m_to_string, scaling_factor):
    """
    @param L: the laplacian corresponding to tips of the tree
    @param eigensplit: the split induced by the fiedler vector
    @param ordered_tip_names: names of the tips of the tree conformant to v and L
    @param m_to_string: a function that converts a matrix to a string
    @param scaling_factor: show the Laplacian scaled by this factor
    @return: a multi-line string
    """
    out = StringIO()
    n = len(L)
    ordered_label_sets = [set([i]) for i in range(n)]
    all_labels = set(range(n))
    for i, child in enumerate(eigensplit):
        complement = all_labels - child
        L_child = SchurAlgebra.mmerge(L, complement) 
        print >> out, 'the Schur complement in the Laplacian of child tree', i+1, 'scaled by', scaling_factor
        print >> out, m_to_string(scaling_factor * L_child)
        print >> out
        child_label_sets = SchurAlgebra.vmerge(ordered_label_sets, complement)
        v_child = BuildTreeTopology.laplacian_to_fiedler(L_child) 
        print >> out, 'the Fiedler split of the Schur complement in the Laplacian of child tree', i+1
        for label_set, value in zip(child_label_sets, v_child):
            s = label_set_to_string(label_set, ordered_tip_names)
            print >> out, s, ':', value
        print >> out
    return out.getvalue().strip()

def get_subtree_messages(D, eigensplit, ordered_tip_names):
    """
    @param D: the matrix of pairwise distances among tips of the tree
    @param eigensplit: the split induced by the fiedler vector
    @param ordered_tip_names: names of the tips of the tree conformant to v and D
    @return: a multi-line string
    """
    out = StringIO()
    n = len(D)
    ordered_label_sets = [set([i]) for i in range(n)]
    all_labels = set(range(n))
    for i, child in enumerate(eigensplit):
        complement = all_labels - child
        D_child = MatrixUtil.get_principal_submatrix(D, list(sorted(child)))
        child_label_sets = SchurAlgebra.vdelete(ordered_label_sets, complement)
        v_child = BuildTreeTopology.edm_to_fiedler(D_child) 
        print >> out, 'the Fiedler split of Schur complements of subtree', i+1
        for label_set, value in zip(child_label_sets, v_child):
            s = label_set_to_string(label_set, ordered_tip_names)
            print >> out, s, ':', value
        print >> out
    return out.getvalue().strip()

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    out = StringIO()
    # get the tree
    tree = NewickIO.parse(fs.tree, FelTree.NewickTree)
    # assert that each node is named
    for node in tree.preorder():
        if not node.name:
            raise HandlingError('each node in the tree must have a name')
    # get the function that converts a matrix to a string
    if fs.plain_matrix:
        m_to_string = MatrixUtil.m_to_string
    elif fs.latex_matrix:
        m_to_string = MatrixUtil.m_to_latex_string
    # print the results for the split of the full tree
    print >> out, get_full_tree_message(tree, m_to_string)
    print >> out
    # get the alphabetically ordered names of the tips
    ordered_tip_names = list(sorted(tip.get_name() for tip in tree.gen_tips()))
    # get the corresponding ordered ids
    tip_name_to_id = dict((tip.get_name(), id(tip)) for tip in tree.gen_tips())
    ordered_tip_ids = [tip_name_to_id[name] for name in ordered_tip_names]
    # get the distance matrix defined by the tips of the tree
    D = np.array(tree.get_partial_distance_matrix(ordered_tip_ids))
    L = Euclid.edm_to_laplacian(D)
    #print >> out, 'the Laplacian obtained from the full tree by Schur complementation:'
    #print >> out, MatrixUtil.m_to_string(L)
    #print >> out
    print >> out, 'the Schur complement in the Laplacian of the full tree scaled by', fs.scaling_factor
    print >> out, m_to_string(fs.scaling_factor * L)
    print >> out
    #L_merged = SchurAlgebra.mmerge(L, set([3,4,5]))
    #print >> out, 'the merged Laplacian:'
    #print >> out, MatrixUtil.m_to_string(L_merged)
    #print >> out
    # get the Fiedler cut of the Schur Laplacian
    v = BuildTreeTopology.laplacian_to_fiedler(L)
    eigensplit = BuildTreeTopology.eigenvector_to_split(v)
    print >> out, 'the Fiedler split of the Schur complement of the full tree:'
    for name, value in zip(ordered_tip_names, v):
        print >> out, name, ':', value
    print >> out
    # get the Fiedler cuts of Schur complements of child trees
    print >> out, get_child_messages(L, eigensplit, ordered_tip_names, m_to_string, fs.scaling_factor)
    print >> out
    # get the Fiedler cuts of Schur complements of subtrees
    print >> out, get_subtree_messages(D, eigensplit, ordered_tip_names)
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
