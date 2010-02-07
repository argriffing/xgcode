"""Find a 4-taxon tree that is in some sense equivalent to a partitioned tree.

Find a 4-taxon tree that is comparable to a tree 
whose taxa are partitioned into four groups.
"""

from StringIO import StringIO

import numpy as np

from SnippetUtil import HandlingError
import SnippetUtil
import Form
import Util
import NewickIO
import FelTree
import Newick
import Clustering
import MatrixUtil

def get_block_structure(taxa_a1, taxa_a2, taxa_b1, taxa_b2):
    a1 = [0]*len(taxa_a1)
    a2 = [1]*len(taxa_a2)
    b1 = [2]*len(taxa_b1)
    b2 = [3]*len(taxa_b2)
    return a1 + a2 + b1 + b2

def get_form():
    """
    @return: a list of form objects
    """
    form_objects = [
            Form.MultiLine('subtree_a', 'rooted subtree A',
                '(a:1, (b:2, c:3):4);'),
            Form.MultiLine('taxa_a1', 'first group of taxa in subtree A',
                '\n'.join('a')),
            Form.MultiLine('taxa_a2', 'second group of taxa in subtree A',
                '\n'.join('bc')),
            Form.MultiLine('subtree_b', 'rooted subtree B',
                '((d:5, e:6):9, (f:7, g:8):10);'),
            Form.MultiLine('taxa_b1', 'first group of taxa in subtree B',
                '\n'.join('de')),
            Form.MultiLine('taxa_b2', 'second group of taxa in subtree B',
                '\n'.join('fg')),
            Form.Float('blen', 'branch length between subtree roots',
                10)]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # read the values from the form
    subtree_a = NewickIO.parse(fs.subtree_a, Newick.NewickTree)
    taxa_a1 = Util.get_stripped_lines(StringIO(fs.taxa_a1))
    taxa_a2 = Util.get_stripped_lines(StringIO(fs.taxa_a2))
    subtree_b = NewickIO.parse(fs.subtree_b, Newick.NewickTree)
    taxa_b1 = Util.get_stripped_lines(StringIO(fs.taxa_b1))
    taxa_b2 = Util.get_stripped_lines(StringIO(fs.taxa_b2))
    connecting_branch_length = fs.blen
    # assert that no group of taxa contains duplicates
    for taxa in (taxa_a1, taxa_a2, taxa_b1, taxa_b2):
        if len(set(taxa)) != len(taxa):
            raise HandlingError('one of the lists of taxa contains duplicates')
    # assert that each subtree has at least two tips and no duplicates
    for tree in (subtree_a, subtree_b):
        tip_names = list(node.get_name() for node in tree.gen_tips())
        if len(tip_names) < 2:
            raise HandlingError('each subtree should have at least two tips')
        if len(set(tip_names)) != len(tip_names):
            raise HandlingError('a subtree has duplicate tip names')
    # assert that the partitions are valid
    first_group = ('A', subtree_a, taxa_a1, taxa_a2) 
    second_group = ('B', subtree_b, taxa_b1, taxa_b2)
    for tree_name, tree, taxa_1, taxa_2 in (first_group, second_group):
        tip_names = set(node.get_name() for node in tree.gen_tips())
        for group_name, taxa in (('1', taxa_1), ('2', taxa_2)):
            nonsense_names = list(set(taxa) - set(tip_names))
            msg_a = 'the following taxa in group %s ' % group_name
            msg_b = 'of subtree %s ' % tree_name
            msg_c = 'are not valid tips: %s' % str(nonsense_names)
            message = msg_a + msg_b + msg_c
            if nonsense_names:
                raise HandlingError(message)
        if set(taxa_1) & set(taxa_2):
            msg_a = 'the taxon lists for subtree %s ' % tree_name
            msg_b = 'are not disjoint'
            raise HandlingError(msg_a + msg_b)
        if set(taxa_1) | set(taxa_2) < tip_names:
            msg_a = 'a tip in subtree %s ' % tree_name
            msg_b = 'is not represented in either of the groups'
            raise HandlingError(msg_a + msg_b)
    # define the response
    out = StringIO()
    # get the results for the first method
    do_first_method(subtree_a, subtree_b, taxa_a1, taxa_a2,
            taxa_b1, taxa_b2, connecting_branch_length, out)
    # define the entire tree by connecting the subtrees
    subtree_b.get_root().set_branch_length(connecting_branch_length)
    subtree_a.get_root().add_child(subtree_b.get_root())
    tree = subtree_a
    # define the order and structure of the distance matrix
    block_structure = get_block_structure(taxa_a1, taxa_a2, taxa_b1, taxa_b2)
    name_order = taxa_a1 + taxa_a2 + taxa_b1 + taxa_b2
    # get the distance matrix
    fel_tree = NewickIO.parse(NewickIO.get_newick_string(tree),
            FelTree.NewickTree)
    D = fel_tree.get_distance_matrix(name_order)
    # get the R matrix
    R = Clustering.get_R_balaji(D)
    # get the sums of block elements of R
    block_R = [[0]*4 for i in range(4)]
    for i, block_i in enumerate(block_structure):
        for j, block_j in enumerate(block_structure):
            block_R[block_i][block_j] += R[i][j]
    # show the results from the second method
    do_second_method(fel_tree, taxa_a1, taxa_a2, taxa_b1, taxa_b2, out)
    # show the results from the third method
    tree_m3_a = NewickIO.parse(fs.subtree_a, Newick.NewickTree)
    tree_m3_b = NewickIO.parse(fs.subtree_b, Newick.NewickTree)
    for t in (tree_m3_a, tree_m3_b):
        neo = Newick.NewickNode()
        neo.name = 'special'
        neo.blen = connecting_branch_length / 2
        t.get_root().add_child(neo)
    feltree_m3_a = NewickIO.parse(NewickIO.get_newick_string(tree_m3_a),
            FelTree.NewickTree)
    feltree_m3_b = NewickIO.parse(NewickIO.get_newick_string(tree_m3_b),
            FelTree.NewickTree)
    tree_m3_a = NewickIO.parse(fs.subtree_a, Newick.NewickTree)
    tree_m3_b = NewickIO.parse(fs.subtree_b, Newick.NewickTree)
    new_root = Newick.NewickNode()
    tree_m3_a.get_root().blen = connecting_branch_length / 2
    tree_m3_b.get_root().blen = connecting_branch_length / 2
    new_root.add_child(tree_m3_a.get_root())
    new_root.add_child(tree_m3_b.get_root())
    tree_m3 = Newick.NewickTree(new_root)
    feltree_m3 = NewickIO.parse(NewickIO.get_newick_string(tree_m3),
            FelTree.NewickTree)
    branch_d2 = connecting_branch_length / 2
    do_third_method(feltree_m3_a, feltree_m3_b, feltree_m3,
            branch_d2, taxa_a1, taxa_a2, taxa_b1, taxa_b2, out)
    # show the expected results
    print >> out, 'M:'
    print >> out, MatrixUtil.m_to_string(R)
    print >> out, 'M summed within blocks:'
    print >> out, MatrixUtil.m_to_string(block_R)
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()

def do_third_method(tree_a, tree_b, tree, branch_d2,
        taxa_a1, taxa_a2, taxa_b1, taxa_b2, out):
    print >> out, 'third method:'
    # get the covariance matrices of the mini trees
    ordered_names = taxa_a1 + taxa_a2 + taxa_b1 + taxa_b2
    ordered_names_a = taxa_a1 + taxa_a2 + ['special']
    ordered_names_b = taxa_b1 + taxa_b2 + ['special']
    block_structure = get_block_structure(taxa_a1, taxa_a2, taxa_b1, taxa_b2)
    block_structure_a = [0]*len(taxa_a1) + [1]*len(taxa_a2) + [2]
    block_structure_b = [0]*len(taxa_b1) + [1]*len(taxa_b2) + [2]
    cov = np.array(tree.get_covariance_matrix(ordered_names))
    cov_a = np.array(tree_a.get_covariance_matrix(ordered_names_a))
    cov_b = np.array(tree_b.get_covariance_matrix(ordered_names_b))
    prec_a = np.linalg.inv(cov_a)
    block_prec_a = [[0]*3 for i in range(3)]
    for i, block_i in enumerate(block_structure_a):
        for j, block_j in enumerate(block_structure_a):
            block_prec_a[block_i][block_j] += prec_a[i][j]
    prec_b = np.linalg.inv(cov_b)
    block_prec_b = [[0]*3 for i in range(3)]
    for i, block_i in enumerate(block_structure_b):
        for j, block_j in enumerate(block_structure_b):
            block_prec_b[block_i][block_j] += prec_b[i][j]
    a = block_prec_a[0][0]
    b = block_prec_a[0][1]
    d = block_prec_a[1][1]
    e = block_prec_b[0][0]
    f = block_prec_b[0][1]
    h = block_prec_b[1][1]
    x = branch_d2
    # make the block M matrix using a clever formula
    Q_a = [
            [a, b, 0, 0],
            [b, d, 0, 0],
            [0, 0, e, f],
            [0, 0, f, h]]
    den_a = (a + 2*b + d + 1/x)
    den_b = (e + 2*f + h + 1/x)
    Q_b = [
            [(a+b)*(a+b)/den_a, (a+b)*(b+d)/den_a, 0, 0],
            [(b+d)*(a+b)/den_a, (b+d)*(b+d)/den_a, 0, 0],
            [0, 0, (e+f)*(e+f)/den_b, (e+f)*(f+h)/den_b],
            [0, 0, (f+h)*(e+f)/den_b, (f+h)*(f+h)/den_b]]
    glom = a+2*b+d+e+2*f+h+2*x*(a+2*b+d)*(e+2*f+h)
    den_a2 = (den_a / den_b) * glom
    den_b2 = (den_b / den_a) * glom
    Q_c = [
            [(a+b)*(a+b)/den_a2, (a+b)*(b+d)/den_a2,
                (a+b)*(e+f)/glom, (a+b)*(f+h)/glom],
            [(b+d)*(a+b)/den_a2, (b+d)*(b+d)/den_a2,
                (b+d)*(e+f)/glom, (b+d)*(f+h)/glom],
            [(e+f)*(a+b)/glom, (e+f)*(b+d)/glom,
                (e+f)*(e+f)/den_b2, (e+f)*(f+h)/den_b2],
            [(f+h)*(a+b)/glom, (f+h)*(b+d)/glom,
                (f+h)*(e+f)/den_b2, (f+h)*(f+h)/den_b2]]
    Q = np.array(Q_a) - np.array(Q_b) - np.array(Q_c)
    print >> out, 'cleverly constructed block M:'
    print >> out, MatrixUtil.m_to_string(Q)
    # make the equivalent tree
    a_star = (b+d)/(a*d-b*b)
    b_star = (a+b)/(a*d-b*b)
    c_star = (f+h)/(e*h-f*f)
    d_star = (e+f)/(e*h-f*f)
    e_star = 2*x - b/(a*d-b*b) - f/(e*h-f*f)
    print >> out, 'using the block precision matrix:'
    print >> out, 'equivalent subtree A:', '(a1:%f, a2:%f);' % (a_star, b_star)
    print >> out, 'equivalent subtree B:', '(b1:%f, b2:%f);' % (c_star, d_star)
    print >> out, 'equivalent connecting branch length:', e_star
    # make the block M matrix using Eric's formula (corrected)
    A, B, C, D, E = a_star, b_star, c_star, d_star, e_star
    H = A*B*(C+D) + C*D*(A+B) + E*(C+D)*(A+B)
    Q = [[1/H]*4 for i in range(4)]
    Q[0][0] *= B*D + B*C + C*D + E*(C+D)
    Q[0][1] *= -C*D - E*(C+D)
    Q[0][2] *= -B*D
    Q[0][3] *= -B*C
    Q[1][0] *= -C*D - E*(C+D)
    Q[1][1] *= A*D + A*C + C*D + E*(C+D)
    Q[1][2] *= -A*D
    Q[1][3] *= -A*C
    Q[2][0] *= -B*D
    Q[2][1] *= -A*D
    Q[2][2] *= A*B + B*D + A*D + E*(A+B)
    Q[2][3] *= -A*B - E*(A+B)
    Q[3][0] *= -B*C
    Q[3][1] *= -A*C
    Q[3][2] *= -A*B - E*(A+B)
    Q[3][3] *= A*B + B*C + A*C + E*(A+B)
    print >> out, 'reconstructed block M:'
    print >> out, MatrixUtil.m_to_string(Q)
    M = Clustering.get_R_balaji(cov)
    M_a = Clustering.get_R_balaji(cov_a)
    M_b = Clustering.get_R_balaji(cov_b)
    block_M = [[0]*4 for i in range(4)]
    for i, block_i in enumerate(block_structure):
        for j, block_j in enumerate(block_structure):
            block_M[block_i][block_j] += M[i][j]
    block_M_a = [[0]*3 for i in range(3)]
    for i, block_i in enumerate(block_structure_a):
        for j, block_j in enumerate(block_structure_a):
            block_M_a[block_i][block_j] += M_a[i][j]
    block_M_b = [[0]*3 for i in range(3)]
    for i, block_i in enumerate(block_structure_b):
        for j, block_j in enumerate(block_structure_b):
            block_M_b[block_i][block_j] += M_b[i][j]
    c_1 = block_M_a[0][1]
    c_2 = block_M_a[0][2]
    c_3 = block_M_a[1][2]
    denominator = (c_1*c_2) + (c_2*c_3) + (c_3*c_1)
    a_star = -c_3 / denominator
    b_star = -c_2 / denominator
    e_star_a = -c_1 / denominator
    c_1 = block_M_b[0][1]
    c_2 = block_M_b[0][2]
    c_3 = block_M_b[1][2]
    denominator = (c_1*c_2) + (c_2*c_3) + (c_3*c_1)
    c_star = -c_3 / denominator
    d_star = -c_2 / denominator
    e_star_b = -c_1 / denominator
    e_star = e_star_a + e_star_b
    print >> out, 'using the block M matrix:'
    print >> out, 'equivalent subtree A:', '(a1:%f, a2:%f);' % (a_star, b_star)
    print >> out, 'equivalent subtree B:', '(b1:%f, b2:%f);' % (c_star, d_star)
    print >> out, 'equivalent connecting branch length:', e_star
    print >> out, 'calculated block M:'
    print >> out, MatrixUtil.m_to_string(block_M)
    print >> out

def do_second_method(tree, taxa_a1, taxa_a2, taxa_b1, taxa_b2, out):
    # get the covariance matrix
    ordered_names = taxa_a1 + taxa_a2 + taxa_b1 + taxa_b2
    cov = np.array(tree.get_covariance_matrix(ordered_names))
    # invert the covariance matrix to make the precision matrix
    prec = np.linalg.inv(cov)
    # take the block sums of the precision matrix
    block_structure = get_block_structure(taxa_a1, taxa_a2, taxa_b1, taxa_b2)
    name_order = taxa_a1 + taxa_a2 + taxa_b1 + taxa_b2
    block_prec = [[0]*4 for i in range(4)]
    for i, block_i in enumerate(block_structure):
        for j, block_j in enumerate(block_structure):
            block_prec[block_i][block_j] += prec[i][j]
    # invert the block summed precision matrix
    reduced_cov = np.linalg.inv(np.array(block_prec))
    # extract the branch lengths from the reduced covariance matrix
    a = reduced_cov[0][0] - reduced_cov[0][1]
    b = reduced_cov[1][1] - reduced_cov[0][1]
    c = reduced_cov[2][2] - reduced_cov[2][3]
    d = reduced_cov[3][3] - reduced_cov[2][3]
    e = reduced_cov[0][1] + reduced_cov[2][3]
    # define the distance matrix for the reduced tree
    reduced_D = [
            [0, a+b, a+e+c, a+e+d],
            [b+a, 0, b+e+c, b+e+d],
            [c+e+a, c+e+b, 0, c+d],
            [d+e+a, d+e+b, d+c, 0]]
    # get the R matrix of the reduced tree
    reduced_R = Clustering.get_R_balaji(reduced_D)
    print >> out, 'second method:'
    print >> out, 'equivalent subtree A:', '(a1:%f, a2:%f);' % (a, b)
    print >> out, 'equivalent subtree B:', '(b1:%f, b2:%f);' % (c, d)
    print >> out, 'equivalent connecting branch length:', e
    print >> out, 'M for the equivalent tree:'
    print >> out, MatrixUtil.m_to_string(reduced_R)
    print >> out

def do_first_method(subtree_a, subtree_b,
        taxa_a1, taxa_a2, taxa_b1, taxa_b2, connecting_branch_length, out):
    # define the branch lengths of the reduced tree
    blen_a1, blen_a2, blen_ar = get_branch_length_equivalents(
            subtree_a, taxa_a1, taxa_a2)
    blen_b1, blen_b2, blen_br = get_branch_length_equivalents(
            subtree_b, taxa_b1, taxa_b2)
    # define the distance matrix of the reduced tree
    a, b = blen_a1, blen_a2
    c, d = blen_b1, blen_b2
    e = connecting_branch_length + blen_ar + blen_br
    reduced_D = [
            [0, a+b, a+e+c, a+e+d],
            [b+a, 0, b+e+c, b+e+d],
            [c+e+a, c+e+b, 0, c+d],
            [d+e+a, d+e+b, d+c, 0]]
    # get the R matrix of the reduced tree
    reduced_R = Clustering.get_R_balaji(reduced_D)
    print >> out, 'first method:'
    print >> out, 'equivalent subtree A:', '(a1:%f, a2:%f);' % (a, b)
    print >> out, 'equivalent subtree B:', '(b1:%f, b2:%f);' % (c, d)
    print >> out, 'equivalent connecting branch length:', e
    print >> out, 'M for the equivalent tree:'
    print >> out, MatrixUtil.m_to_string(reduced_R)
    print >> out

def get_branch_length_equivalents(tree, first_taxa, second_taxa):
    """
    @param tree: a newick tree
    @param first_taxa: a set of tip names
    @param second_taxa: another set of tip names
    @return: a triple (first distance, second distance, root distance)
    """
    # get the root-augmented distance matrices
    D_aug = get_root_augmented_distance_matrix(tree, first_taxa, second_taxa)
    # get the R matrix
    R_aug = Clustering.get_R_balaji(D_aug)
    # Get the matrix whose elements are block element sums
    # of the root-augmented R matrix.
    block_structure = [0]*len(first_taxa) + [1]*len(second_taxa) + [2]
    B = [[0]*3 for i in range(3)]
    for i, block_i in enumerate(block_structure):
        for j, block_j in enumerate(block_structure):
            B[block_i][block_j] += R_aug[i][j]
    # get the new branch lengths for the subtree
    denominator = 2 * (B[0][1]*B[1][2] + B[1][2]*B[2][0] + B[2][0]*B[0][1])
    blen_first = B[1][2] / denominator
    blen_second = B[2][0] / denominator
    blen_root = B[0][1] / denominator
    return blen_first, blen_second, blen_root

def get_root_augmented_distance_matrix(tree_in, first_taxa, second_taxa):
    """
    @param tree_in: a newick tree
    @param first_taxa: a set of tip names
    @param second_taxa: another set of tip names
    @return: a distance matrix
    """
    # first convert the tree to the appropriate data structure
    tree = NewickIO.parse(NewickIO.get_newick_string(tree_in),
            FelTree.NewickTree)
    # now get the ordered ids
    ordered_ids = []
    for taxa in (first_taxa, second_taxa):
        for node in tree.gen_tips():
            if node.get_name() in taxa:
                ordered_ids.append(id(node))
    ordered_ids.append(id(tree.get_root()))
    # now get the distance matrix
    return tree.get_partial_distance_matrix(ordered_ids)

