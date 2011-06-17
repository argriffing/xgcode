"""
Show leaves in nodal partitions of a tree, for full and leaf-only matrices.
"""

from StringIO import StringIO

import numpy as np
import scipy

import Form
import FormOut
import Ftree
import FtreeIO
import MatrixUtil
import Newick
import forest
import const

g_tree_string = const.read('20110616a').rstrip()


def get_form():
    """
    @return: a list of form objects
    """
    # reformat the tree string
    tree = Newick.parse(g_tree_string, Newick.NewickTree)
    formatted_tree_string = Newick.get_narrow_newick_string(tree, 75)
    # define the form objects
    form_objects = [
            Form.MultiLine('tree_string', 'newick tree',
                formatted_tree_string)]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_harmonically_extended_MDS(T, B, leaves, internal):
    """
    Use harmonically extended 2D MDS.
    """
    Lbb = Ftree.TB_to_L_block(T, B, internal, internal)
    Lba = Ftree.TB_to_L_block(T, B, internal, leaves)
    L_schur = Ftree.TB_to_L_schur(T, B, leaves)
    W, V = scipy.linalg.eigh(L_schur, eigvals=(1, 2))
    V = V * np.reciprocal(np.sqrt(W))
    Y = -np.dot(np.dot(np.linalg.pinv(Lbb), Lba), V)
    return np.vstack([V, Y])

def get_full_distance_MDS(T, B, vertices):
    D_full = Ftree.TB_to_D(T, B, vertices)
    G_neg = MatrixUtil.double_centered(D_full) / 2.0
    W_neg, V = scipy.linalg.eigh(G_neg, eigvals=(0,1))
    return V * np.sqrt(np.abs(W_neg))

def get_fp1_ordered_leaf_partition(T, v_to_value):
    """
    The notation fp1 means Fiedler plus one.
    Note that a sequence is returned, so the partition is ordered.
    @return: a sequence of three leaf sets
    """
    nvertices = len(v_to_value)
    # find the two edges where the valuation changes at the endpoints
    cuts = [frozenset([a,b]) for a, b in T if v_to_value[a]*v_to_value[b]<0]
    if len(cuts) != 2:
        msg_a = 'expected exactly two nodes (in the sense of Courant) '
        msg_b = 'but found ' + str(len(cuts))
        raise ValueError(msg_a + msg_b)
    # create a forest by cauterizing the cut edges
    next_vertex = nvertices
    T_cut = set(T)
    for cut in cuts:
        a, b = cut
        T_cut.remove(cut)
        T_cut.add(frozenset([a, next_vertex]))
        T_cut.add(frozenset([b, next_vertex+1]))
        next_vertex += 2
    # get the unordered vertex partition from the forest
    partition = list(forest.VT_to_vertex_partition(set(v_to_value), T_cut))
    # count the number of cauterization stubs per vertex set
    nstubs_list = [
            sum(1 for v in p if v >= nvertices) for p in partition]
    # reorder the partition according to the number of cauterization stubs
    pairs = sorted(zip(nstubs_list, partition))
    pairs[0], pairs[1] = pairs[1], pairs[0]
    partition = [v_set for nstubs, v_set in pairs]
    return partition

def get_2D_report(N, leaf_set, title, fiedler_partition, fp1_partition):
    out = StringIO()
    print >> out, title
    print >> out, 'Fiedler nodal domain leaf partition:'
    for prefix, v_set in zip(('A', 'B'), fiedler_partition):
        for v in sorted(v_set & leaf_set):
            print >> out, prefix, N.get(v, '?')
    print >> out
    print >> out, title
    print >> out, 'Fiedler+1 nodal domain leaf partition:'
    for prefix, v_set in zip(('A', 'B', 'C'), fp1_partition):
        for v in sorted(v_set & leaf_set):
            print >> out, prefix, N.get(v, '?')
    return out.getvalue()

def get_response_content(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: the response
    """
    T, B, N = FtreeIO.newick_to_TBN(fs.tree_string)
    leaves = Ftree.T_to_leaves(T)
    internal = Ftree.T_to_internal_vertices(T)
    vertices = leaves + internal
    # compute the 2D MDS of the full tree
    MDS_full = get_full_distance_MDS(T, B, vertices)
    # compute the harmonic extension of 2D MDS of the leaves
    MDS_harmonic = get_harmonically_extended_MDS(T, B, leaves, internal)
    # get the Fiedler MDS leaf partition for full 2D MDS
    v_to_value = dict(zip(vertices, MDS_full[:,0]))
    neg_leaves = frozenset([v for v in leaves if v_to_value[v] < 0])
    pos_leaves = frozenset([v for v in leaves if v_to_value[v] >= 0])
    full_mds_fiedler_partition = [neg_leaves, pos_leaves]
    # get the Fiedler MDS leaf partition for harmonic 2D MDS
    v_to_value = dict(zip(vertices, MDS_harmonic[:,0]))
    neg_leaves = frozenset([v for v in leaves if v_to_value[v] < 0])
    pos_leaves = frozenset([v for v in leaves if v_to_value[v] >= 0])
    harmonic_mds_fiedler_partition = [neg_leaves, pos_leaves]
    # get the Fiedler plus one MDS leaf partition for full 2D MDS
    v_to_value = dict(zip(vertices, MDS_full[:,1]))
    full_mds_fp1_partition = get_fp1_ordered_leaf_partition(T, v_to_value)
    # get the Fiedler plus one MDS leaf partition for harmonic 2D MDS
    v_to_value = dict(zip(vertices, MDS_harmonic[:,1]))
    harmonic_mds_fp1_partition = get_fp1_ordered_leaf_partition(T, v_to_value)
    # write the output
    out = StringIO()
    print >> out, get_2D_report(
            N, set(leaves), 'Full distance 2D MDS',
            full_mds_fiedler_partition, full_mds_fp1_partition)
    print >> out
    print >> out, get_2D_report(
            N, set(leaves), 'Harmonically extended 2D MDS',
            harmonic_mds_fiedler_partition, harmonic_mds_fp1_partition)
    return out.getvalue()

