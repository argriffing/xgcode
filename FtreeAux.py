"""
Auxiliary functions for ftrees.

This is a more functional programming version of SpatialTree layout.
It is less stateful so it might be more flexible.
"""

import unittest
import math
from collections import defaultdict

import numpy as np

import day
import Ftree
import FtreeIO
import LeafWeights

NEG_EDGE = -1
NUL_EDGE = 0
POS_EDGE = 1
ALT_EDGE = 2


def ndot(*args):
    M = args[0]
    for B in args[1:]:
        M = np.dot(M, B)
    return M

#####################################################
# This section is for EqualDaylight layout.


def _build_dtree(dtree, v, v_to_sinks, v_to_location, v_to_dtree_id, count):
    """
    @param dtree: an object that lives in the C extension
    @param v: vertex
    @param v_to_sinks: computed from a rooted Ftree
    @parma v_to_location: vertex to location map
    @param v_to_dtree_id: this is filled as the dtree is built
    @param count: the number of nodes added so far
    """
    x, y = v_to_location[v]
    v_to_dtree_id[v] = count
    count += 1
    dtree.begin_node(v_to_dtree_id[v])
    dtree.set_x(x)
    dtree.set_y(y)
    sinks = v_to_sinks.get(v, [])
    for sink in sorted(sinks):
        count = _build_dtree(
                dtree, sink, v_to_sinks, v_to_location, v_to_dtree_id, count)
    dtree.end_node()
    return count

def equal_daylight_layout(T, B, iteration_count):
    """
    @param T: topology
    @param B: branch lengths
    """
    R = Ftree.T_to_R_canonical(T)
    r = Ftree.R_to_root(R)
    # create the initial equal arc layout
    v_to_location = equal_arc_layout(T, B)
    # use sax-like events to create a parallel tree in the C extension
    v_to_sinks = Ftree.R_to_v_to_sinks(R)
    v_to_dtree_id = {}
    dtree = day.Day()
    count = _build_dtree(
            dtree, r, v_to_sinks, v_to_location, v_to_dtree_id, 0)
    # repeatedly reroot and equalize
    v_to_neighbors = Ftree.T_to_v_to_neighbors(T)
    for i in range(iteration_count):
        for v in Ftree.T_to_inside_out(T):
            neighbor_count = len(v_to_neighbors[v])
            if neighbor_count > 2:
                dtree.select_node(v_to_dtree_id[v])
                dtree.reroot()
                dtree.equalize()
    # extract the x and y coordinates from the dtree
    v_to_location = {}
    for v, dtree_id in v_to_dtree_id.items():
        dtree.select_node(dtree_id)
        x = dtree.get_x()
        y = dtree.get_y()
        v_to_location[v] = (x, y)
    return v_to_location


#####################################################
# This section is for EqualArc layout.


def equal_arc_layout(T, B):
    """
    @param T: tree topology
    @param B: branch lengths
    @return: a map from vertex to location
    """
    # arbitrarily root the tree
    R = Ftree.T_to_R_canonical(T)
    r = Ftree.R_to_root(R)
    # map vertices to subtree tip count
    v_to_sinks = Ftree.R_to_v_to_sinks(R)
    v_to_count = {}
    for v in Ftree.R_to_postorder(R):
        sinks = v_to_sinks.get(v, [])
        if sinks:
            v_to_count[v] = sum(v_to_count[sink] for sink in sinks)
        else:
            v_to_count[v] = 1
    # create the equal arc angles
    v_to_theta = {}
    _force_equal_arcs(
            v_to_sinks, v_to_count, v_to_theta,
            r, -math.pi, math.pi)
    # convert angles to coordinates
    v_to_source = Ftree.R_to_v_to_source(R)
    v_to_location = {}
    _update_locations(
            R, B,
            v_to_source, v_to_sinks, v_to_theta, v_to_location,
            r, (0, 0), 0)
    return v_to_location

def _force_equal_arcs(
        v_to_sinks, v_to_count, v_to_theta,
        v, min_theta, max_theta):
    """
    Use the equal angle method.
    Lay out the tree with non-intersecting branches.
    Define the angles in this subtree to be within the specified range.
    Fill v_to_theta.
    """
    v_to_theta[v] = (min_theta + max_theta) / 2
    sinks = v_to_sinks.get(v, [])
    if sinks:
        subtree_tip_count = v_to_count[v]
        cumulative_theta = 0
        for child in sinks:
            sub_subtree_tip_count = v_to_count[child]
            dtheta = max_theta - min_theta
            aliquot = dtheta * sub_subtree_tip_count / float(subtree_tip_count)
            low = min_theta + cumulative_theta - v_to_theta[v]
            cumulative_theta += aliquot
            high = min_theta + cumulative_theta - v_to_theta[v]
            _force_equal_arcs(
                    v_to_sinks, v_to_count, v_to_theta,
                    child, low, high)

def _update_locations(
        R, B,
        v_to_source, v_to_sinks, v_to_theta, v_to_location,
        v, last_location, last_theta):
    """
    Fill v_to_location.
    @param last_location: an (x, y) pair or None if root
    @param last_theta: the direction we were last facing
    """
    theta = last_theta + v_to_theta[v]
    parent = v_to_source.get(v, None)
    if parent is None:
        branch_length = None
    else:
        u_edge = frozenset((v, parent))
        branch_length = B.get(u_edge, None)
    if branch_length is None:
        v_to_location[v] = last_location
    else:
        x, y = last_location
        dx = branch_length*math.cos(theta)
        dy = branch_length*math.sin(theta)
        v_to_location[v] = (x + dx, y + dy)
    sinks = v_to_sinks.get(v, [])
    for child in sinks:
        _update_locations(R, B,
                v_to_source, v_to_sinks, v_to_theta, v_to_location,
                child, v_to_location[v], theta)


#####################################################
# This section is for leaf weight calculation.

def get_internal_vertex_to_leaf_distn(T, B):
    """
    Return a map from an internal vertex to a leaf distribution.
    @return: a dictionary that maps an internal vertex to a leaf distribution
    """
    leaves = Ftree.T_to_leaves(T)
    internal = Ftree.T_to_internal_vertices(T)
    vertices = leaves + internal
    # Get pieces of the Laplacian matrix.
    Lbb = Ftree.TB_to_L_block(T, B, internal, internal)
    Lba = Ftree.TB_to_L_block(T, B, internal, leaves)
    # Get the numpy array of harmonic extensions to previously unknown values.
    interpolator = -np.dot(np.linalg.pinv(Lbb), Lba)
    #print 'L interpolator:'
    #print interpolator.shape
    #print interpolator
    d = {}
    for i, v in enumerate(internal):
        distn = {}
        for j, leaf in enumerate(leaves):
            distn[leaf] = interpolator[i, j]
        d[v] = distn
    return d

def get_internal_vertex_to_leaf_distn_cov(T, B):
    """
    This is a possibly equivalent formualtion.
    It is based on Schur complementation in the unrooted covariance matrix.
    Return a map from an internal vertex to a leaf distribution.
    @return: a dictionary that maps an internal vertex to a leaf distribution
    """
    leaves = Ftree.T_to_leaves(T)
    internal = Ftree.T_to_internal_vertices(T)
    vertices = leaves + internal
    # Get the full tree Laplacian matrix.
    L = Ftree.TB_to_L_principal(T, B, vertices)
    # Get the unrooted covariance matrix.
    HSH = np.linalg.pinv(L)
    # Use the multivariate normal distribution wikipedia page
    # for conditional distributions.
    nleaves = len(leaves)
    ninternal = len(internal)
    #
    # This interpolator works.
    #Lbb = L[nleaves:, nleaves:]
    #Lba = L[nleaves:, :nleaves]
    #interpolator = -ndot(np.linalg.pinv(Lbb), Lba)
    #
    # This interpolator seems like it should work but it does not.
    Saa = HSH[:nleaves, :nleaves]
    Sba = HSH[nleaves:, :nleaves]
    #print 'det(Saa)'
    #print np.linalg.det(Saa)
    interpolator = ndot(Sba, np.linalg.pinv(Saa))
    #
    # Try a hack.
    #eps = 1e-12
    #nvertices = len(vertices)
    #J = np.ones((nvertices, nvertices))
    #Saa = (HSH + J)[:nleaves, :nleaves]
    #Sba = (HSH + J)[nleaves:, :nleaves]
    #interpolator = ndot(Sba, np.linalg.pinv(Saa))
    #
    #print 'cov interpolator:'
    #print interpolator.shape
    #print interpolator
    d = {}
    for i, v in enumerate(internal):
        distn = {}
        for j, leaf in enumerate(leaves):
            distn[leaf] = interpolator[i, j]
        d[v] = distn
    return d

def get_leaf_distn_acl(R, B):
    """
    This is a possibly equivalent formulation.
    It is based on Felsenstein weights.
    """
    # Get the vertex order.
    T = Ftree.R_to_T(R)
    r = Ftree.R_to_root(R)
    leaves = Ftree.T_to_leaves(T)
    non_r_internal = [v for v in Ftree.T_to_internal_vertices(T) if v != r]
    vertices = leaves + non_r_internal + [r]
    # Get the pseudoinverse of the Laplacian.
    # This is also the doubly centered covariance matrix.
    L = Ftree.TB_to_L_principal(T, B, vertices)
    HSH = np.linalg.pinv(L)
    # Decenter the covariance matrix using the root.
    # This should give the rooted covariance matrix
    # which is M in the appendix of Weights for Data Related by a Tree
    # by Altschul, Carroll, and Lipman, 1989.
    e = np.ones_like(HSH[-1])
    J = np.ones_like(HSH)
    M = HSH - np.outer(e, HSH[-1]) - np.outer(HSH[-1], e) + HSH[-1,-1]*J
    # Pick out the part corresponding to leaves.
    nleaves = len(leaves)
    S = M[:nleaves, :nleaves]
    S_pinv = np.linalg.pinv(S)
    # Normalized row or column sums of inverse of M gives the leaf distribution.
    w = S_pinv.sum(axis=0) / S_pinv.sum()
    return dict((v, w[i]) for i, v in enumerate(leaves))

def get_leaf_distn_schur(R, B):
    """
    This is a possibly equivalent formulation.
    It is based on removing all internal vertices except the root
    by Schur complement.
    """
    # Get the vertex order.
    # This order is different from the acl order.
    T = Ftree.R_to_T(R)
    r = Ftree.R_to_root(R)
    leaves = Ftree.T_to_leaves(T)
    non_r_internal = [v for v in Ftree.T_to_internal_vertices(T) if v != r]
    vertices = leaves + [r] + non_r_internal
    # Get the combinatorial Laplacian matrix
    # and Schur complement out all of the non-root internal vertices.
    L_schur = Ftree.TB_to_L_schur(T, B, leaves + [r])
    # Get the vector of negative weights between the root and the leaves.
    w_unnormalized = L_schur[-1, :-1]
    # Get the normalized weight vector
    w = w_unnormalized / w_unnormalized.sum()
    return dict((v, w[i]) for i, v in enumerate(leaves))


#####################################################
# This section is for splitting the branches.


def values_to_color(value_a, value_b, eps):
    """
    The color is really an integer related to a combination of vertex signs.
    Edges with color 2 should be split before visualization.
    Edges with color -1 are drawn thin.
    Edges with color 0 are drawn wavy.
    Edges with color 1 are drawn thick.
    Color -1: one vertex is very negative and the other is not very positive.
    Color 0: both vertices have negligible valuation.
    Color 1: one vertex is very positive and the other is not very negative.
    Color 2: one vertex is very positive and the other is very negative.
    @param value_a: value at one vertex
    @param value_b: value at the other vertex
    @param eps: small positive
    @return: a color in -1, 0, 1, 2
    """
    low, high = sorted((value_a, value_b))
    if low < -eps and high < eps:
        return NEG_EDGE
    elif high > eps and low > -eps:
        return POS_EDGE
    elif low > -eps and high < eps:
        return NUL_EDGE
    else:
        return ALT_EDGE

def break_branches_by_vertex_sign(T, B, v_to_value, eps):
    """
    This function modifies T, B, and v_to_value.
    I guess that means it is not very functional programming.
    Add degree two vertices.
    After the degree two vertices have been added,
    each branch should satisfy exactly one of the following three conditions.
    Both endpoints are near zero.
    One endpoint is very positive and the other is not very negative.
    One endpoint is very negative and the other is not very positive.
    @param T: topology
    @param B: branch lengths
    @parm v_to_value: map a vertex to a valuation
    @param eps: small
    """
    vertices = Ftree.T_to_order(T)
    next_v = max(vertices) + 1
    d_edges = Ftree.T_to_outside_in_edges(T)
    for d_edge in d_edges:
        va, vb = d_edge
        vsrc = v_to_value[va]
        vdst = v_to_value[vb]
        if values_to_color(vsrc, vdst, eps) != ALT_EDGE:
            continue
        # find the crossing point
        t = -vsrc / (vdst - vsrc)
        # get the branch length
        u_edge = frozenset(d_edge)
        branch_length_a_b = B[u_edge]
        # remove the old edge
        T.remove(u_edge)
        del B[u_edge]
        # add the new edges
        u_a_mid = frozenset((va, next_v))
        u_b_mid = frozenset((vb, next_v))
        T.add(u_a_mid)
        T.add(u_b_mid)
        B[u_a_mid] = branch_length_a_b * t
        B[u_b_mid] = branch_length_a_b * (1.0 - t)
        # update the value for the new vertex
        v_to_value[next_v] = 0.0
        # increment the vertex counter
        next_v += 1


#####################################################
# The remaining functions are miscellaneous.


def harmonically_interpolate(T, B, v_to_value):
    """
    Use the harmonic extension to augment the v_to_value map.
    The T and B data is not modified.
    """
    vertices = Ftree.T_to_order(T)
    # Define the lists of vertices for which the values are known and unknown.
    known = sorted(v_to_value)
    unknown = sorted(set(vertices) - set(v_to_value))
    # If everything is known then we do not need to interpolate.
    if not unknown:
        return
    # Get pieces of the Laplacian matrix.
    Lbb = Ftree.TB_to_L_block(T, B, unknown, unknown)
    Lba = Ftree.TB_to_L_block(T, B, unknown, known)
    # Get the numpy array of known values.
    v_known = np.array([v_to_value[v] for v in known])
    # Get the numpy array of harmonic extensions to previously unknown values.
    v_unknown = -np.dot(np.dot(np.linalg.pinv(Lbb), Lba), v_known)
    # Put the interpolated values into the dictionary.
    for vertex, value in zip(unknown, v_unknown):
        v_to_value[vertex] = value

def color_edges_by_vertex_sign(T, v_to_value, eps):
    """
    Colors should be integers -1, 0, 1.
    @param T: topology
    @param v_to_value: vertex to value
    @param eps: small positive
    @return: edge_to_color
    """
    edge_to_color = {}
    for u_edge in Ftree.T_to_edges(T):
        va, vb = u_edge
        a_val, b_val = v_to_value[va], v_to_value[vb]
        edge_to_color[u_edge] = values_to_color(a_val, b_val, eps)
    return edge_to_color

def T_to_edge_to_neighbor_edges(T):
    """
    @param T: topology
    @return: map from undirected edge to list of undirected edges
    """
    edge_to_neighbor_edges = defaultdict(list)
    v_to_neighbors = Ftree.T_to_v_to_neighbors(T)
    for edge in Ftree.T_to_edges(T):
        for source in edge:
            for sink in v_to_neighbors[source]:
                if sink not in edge:
                    n_edge = frozenset((source, sink))
                    edge_to_neighbor_edges[edge].append(n_edge)
    return edge_to_neighbor_edges

def get_multi_edges(T, edge_to_color):
    """
    Here a multi-edge is a tuple of vertices, and it is basically a path.
    Each sequential pair of vertices in the multi-edge is an edge of T,
    and each sequential triple of vertices represents two connected edges.
    Each edge represented in a multi-edge must have the same edge color.
    The idea is to concatenate some edges into monochromatic paths.
    This function partitions the edges into monochromatic paths.
    The partition is somewhat arbitrary.
    @param T: topology
    @param edge_to_color: maps an edge to a color
    @return: vertex tuples
    """
    multi_edges = []
    visited_edges = set()
    v_to_neighbors = Ftree.T_to_v_to_neighbors(T)
    edge_to_neighbor_edges = T_to_edge_to_neighbor_edges(T)
    for u_edge in Ftree.T_to_edges(T):
        if u_edge in visited_edges:
            continue
        c = edge_to_color[u_edge]
        # find a neighbor edge of the same color
        target_edge = None
        for n_edge in edge_to_neighbor_edges[u_edge]:
            if n_edge in visited_edges:
                continue
            if edge_to_color[n_edge] == c:
                target_edge = n_edge
                break
        if target_edge:
            # if a target edge was found then add the three vertex path
            vm, = u_edge & n_edge
            va, vb = (u_edge | n_edge) - (u_edge & n_edge)
            edge_amb = (va, vm, vb)
            edge_bma = (vb, vm, va)
            directed_multiedge = min(edge_amb, edge_bma)
            multi_edges.append(directed_multiedge)
            visited_edges.add(n_edge)
            visited_edges.add(u_edge)
        else:
            # if a target edge was not found then add the two vertex path
            directed_multiedge = tuple(sorted(u_edge))
            multi_edges.append(directed_multiedge)
            visited_edges.add(u_edge)
    return multi_edges


class TestFtreeAux(unittest.TestCase):

    def setUp(self):
        np.set_printoptions(linewidth=300)

    def test_get_multi_edges(self):
        T = set(frozenset(pair) for pair in [
            (0,2), (1,2), (2,3), (3,4), (3,5), (5,6), (5,7)])
        edge_to_color = dict((frozenset(pair), c) for pair, c in [
            ((0,2), 'a'), ((1,2), 'a'),
            ((2,3), 'b'),
            ((3,4), 'a'), ((3,5), 'a'),
            ((5,6), 'b'), ((5,7), 'b')])
        observed = set(get_multi_edges(T, edge_to_color))
        expected = set([
            (0,2,1), (2,3), (4,3,5), (6,5,7)])
        self.assertEqual(observed, expected)

    def test_break_branches_by_vertex_sign(self):
        eps = 1e-8
        T = set(frozenset(pair) for pair in [
            (0,1), (0,2), (0,3)])
        B = dict((frozenset(pair), float(b)) for pair, b in [
            ((0,1), 1), ((0,2), 1), ((0,3), 1)])
        v_to_value = {
                0:1.0, 1:1.0, 2:1.0, 3:-3.0}
        break_branches_by_vertex_sign(T, B, v_to_value, eps)
        T_expected = set(frozenset(pair) for pair in [
            (0,1), (0,2), (0,4), (4,3)])
        B_expected = dict((frozenset(pair), float(b)) for pair, b in [
            ((0,1), 1), ((0,2), 1), ((0,4), 0.25), ((4,3), 0.75)])
        v_to_value_expected = {
                0:1.0, 1:1.0, 2:1.0, 3:-3.0, 4:0.0}
        self.assertEqual(T, T_expected)
        self.assertEqual(B, B_expected)
        self.assertEqual(v_to_value, v_to_value_expected)

    def test_harmonically_interpolate(self):
        T = set(frozenset(pair) for pair in [
            (0,1), (0,2), (0,4), (4,3)])
        B = dict((frozenset(pair), float(b)) for pair, b in [
            ((0,1), 1), ((0,2), 1), ((0,4), 0.25), ((4,3), 0.75)])
        v_to_x_offset = {
                0:5.0, 1:7.0, 2:8.0, 3:9.0}
        harmonically_interpolate(T, B, v_to_x_offset)
        v_to_x_offset_expected = {
                0:5.0, 1:7.0, 2:8.0, 3:9.0, 4:6.0}
        self.assertEqual(v_to_x_offset, v_to_x_offset_expected)

    def test_color_edges_by_vertex_sign(self):
        eps = 1e-8
        T = set(frozenset(pair) for pair in [
            (1,2), (2,4), (3,4), (4,5), (5,6), (5,7), (7,8)])
        v_to_value = {
                1:2.0, 2:1.0, 3:1.0, 4:1e-10, 5:0.0, 6:-1.0, 7:-1.0, 8:3.0}
        edge_to_color_observed = color_edges_by_vertex_sign(T, v_to_value, eps)
        edge_to_color_expected = dict((frozenset(p), float(x)) for p, x in [
            ((1,2),POS_EDGE), ((2,4),POS_EDGE), ((3,4),POS_EDGE),
            ((4,5),NUL_EDGE),
            ((5,6),NEG_EDGE), ((5,7),NEG_EDGE), ((7,8),ALT_EDGE)])
        self.assertEqual(edge_to_color_observed, edge_to_color_expected)

    def test_leaf_distn_a(self):
        # Read the example tree.
        example_tree = '(a:2, (b:1, c:1, d:1, e:1)x:1)y;'
        R, B, N = FtreeIO.newick_to_RBN(example_tree)
        T = Ftree.R_to_T(R)
        r = Ftree.R_to_root(R)
        # Get the leaf distribution associated with the root.
        internal_to_leaf_distn = get_internal_vertex_to_leaf_distn(T, B)
        r_to_leaf_distn = internal_to_leaf_distn[r]
        leaves = Ftree.T_to_leaves(T)
        observed_name_weight_pairs = [
                (N[v], r_to_leaf_distn[v]) for v in leaves]
        # Set up the expectation for the test.
        n = 5.0
        expected_name_weight_pairs = []
        expected_first_value = n / (3*n - 2)
        expected_non_first_value = 2 / (3*n - 2)
        expected_name_weight_pairs.append(('a', expected_first_value))
        for name in list('bcde'):
            expected_name_weight_pairs.append((name, expected_non_first_value))
        # Do the comparison for testing.
        expected_d = dict(expected_name_weight_pairs)
        observed_d = dict(observed_name_weight_pairs)
        for v in leaves:
            name = N[v]
            expected_value = expected_d[name]
            observed_value = observed_d[name]
            self.assertTrue(np.allclose(expected_value, observed_value))

    def test_leaf_distn_b(self):
        # Read the example tree.
        example_tree = LeafWeights.g_acl_tree
        R, B, N = FtreeIO.newick_to_RBN(example_tree)
        T = Ftree.R_to_T(R)
        r = Ftree.R_to_root(R)
        # Get the leaf distribution associated with the root.
        internal_to_leaf_distn = get_internal_vertex_to_leaf_distn(T, B)
        leaf_distn = internal_to_leaf_distn[r]
        leaves = Ftree.T_to_leaves(T)
        observed_name_weight_pairs = [
                (N[v], leaf_distn[v]) for v in leaves]
        # Do the comparison for testing.
        observed_name_to_weight = dict(observed_name_weight_pairs)
        for name in LeafWeights.g_acl_ordered_names:
            s_expected = LeafWeights.g_acl_expected_weights[name]
            s_observed = '%.3f' % observed_name_to_weight[name]
            self.assertEqual(s_expected, s_observed)

    def test_leaf_distn_cov(self):
        """
        This test is known to fail.
        """
        # Read the example tree.
        example_tree = LeafWeights.g_acl_tree
        R, B, N = FtreeIO.newick_to_RBN(example_tree)
        T = Ftree.R_to_T(R)
        r = Ftree.R_to_root(R)
        # Get the leaf distribution associated with the root.
        internal_to_leaf_distn = get_internal_vertex_to_leaf_distn_cov(T, B)
        leaf_distn = internal_to_leaf_distn[r]
        leaves = Ftree.T_to_leaves(T)
        observed_name_weight_pairs = [
                (N[v], leaf_distn[v]) for v in leaves]
        # Do the comparison for testing.
        observed_name_to_weight = dict(observed_name_weight_pairs)
        for name in LeafWeights.g_acl_ordered_names:
            s_expected = LeafWeights.g_acl_expected_weights[name]
            s_observed = '%.3f' % observed_name_to_weight[name]
            self.assertEqual(s_expected, s_observed)

    def test_leaf_distn_acl(self):
        # Read the example tree.
        example_tree = LeafWeights.g_acl_tree
        R, B, N = FtreeIO.newick_to_RBN(example_tree)
        T = Ftree.R_to_T(R)
        r = Ftree.R_to_root(R)
        # Get the leaf distribution associated with the root.
        leaf_distn = get_leaf_distn_acl(R, B)
        leaves = Ftree.T_to_leaves(T)
        observed_name_weight_pairs = [
                (N[v], leaf_distn[v]) for v in leaves]
        # Do the comparison for testing.
        observed_name_to_weight = dict(observed_name_weight_pairs)
        for name in LeafWeights.g_acl_ordered_names:
            s_expected = LeafWeights.g_acl_expected_weights[name]
            s_observed = '%.3f' % observed_name_to_weight[name]
            self.assertEqual(s_expected, s_observed)

    def test_leaf_distn_schur(self):
        # Read the example tree.
        example_tree = LeafWeights.g_acl_tree
        R, B, N = FtreeIO.newick_to_RBN(example_tree)
        T = Ftree.R_to_T(R)
        r = Ftree.R_to_root(R)
        # Get the leaf distribution associated with the root.
        leaf_distn = get_leaf_distn_schur(R, B)
        leaves = Ftree.T_to_leaves(T)
        observed_name_weight_pairs = [
                (N[v], leaf_distn[v]) for v in leaves]
        # Do the comparison for testing.
        observed_name_to_weight = dict(observed_name_weight_pairs)
        for name in LeafWeights.g_acl_ordered_names:
            s_expected = LeafWeights.g_acl_expected_weights[name]
            s_observed = '%.3f' % observed_name_to_weight[name]
            self.assertEqual(s_expected, s_observed)


if __name__ == '__main__':
    unittest.main()

