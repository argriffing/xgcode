"""Project the internal tree vertices onto the 2d MDS plane of the leaves.

Project the internal tree vertices onto the 2d plane
defined by the MDS of the leaves.
"""

from StringIO import StringIO
import math
import random

import numpy as np
import cairo

from SnippetUtil import HandlingError
import SnippetUtil
import Form
import FormOut
import CairoUtil
import MatrixUtil
import Clustering
import NewickIO
import FelTree
import Euclid
import TreeSampler
import TreeComparison


#g_tree_string = '(((a:0.05, b:0.05):0.15, c:0.2):0.8, x:1.0, (((m:0.05, n:0.05):0.15, p:0.2):0.8, y:1.0):1.0);'
g_tree_string = '((((a1:2, a2:3):2, a3:1):1, a4:5):4, (((b1:8, b2:2):3, b3:4):1, b4:1):2);'


def get_form():
    """
    @return: a list of form objects
    """
    # define the default tree string
    tree = NewickIO.parse(g_tree_string, FelTree.NewickTree)
    formatted_default_tree_string = NewickIO.get_narrow_newick_string(tree, 60)
    # define the list of form objects
    form_objects = [
            Form.MultiLine('tree', 'tree', formatted_default_tree_string),
            Form.ImageFormat(),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.Image('tree')

def get_rescaled_vector(v):
    """
    @param v: an array or list of floating point values
    @return: a list of floating point values rescaled to be in the range (0, 1)
    """
    width = max(v) - min(v)
    if not width:
        return [.5 for x in v]
    return [(x-min(v)) / width for x in v]

def get_image(projected_points, nleaves, incidence_matrix, ordered_names, width_and_height, image_format):
    """
    @param projected_points: points projected onto the 2d plane
    @param nleaves: the first few points are leaves
    @param incidence_matrix: for drawing connections
    @param ordered_names: the labels corresponding to rows of the row major matrix
    @param width_and_height: the dimensions of the output image
    @param image_format: like 'svg', 'png', 'ps', 'pdf', et cetera
    @return: a string containing the image data
    """
    width, height = width_and_height
    n = len(projected_points)
    # set some values that used to be arguments
    draw_axes = True
    draw_connections = True
    # get eigenvectors scaled to [0, 1]
    va, vb = projected_points.T.tolist()
    rescaled_a = get_rescaled_vector(va)
    rescaled_b = get_rescaled_vector(vb)
    # create the surface
    cairo_helper = CairoUtil.CairoHelper(image_format)
    surface = cairo_helper.create_surface(width, height)
    context = cairo.Context(surface)
    # draw the background
    context.save()
    context.set_source_rgb(.9, .9, .9)
    context.paint()
    context.restore()
    # define the border
    border_fraction = .1
    # draw the axes if requested
    if draw_axes:
        # begin drawing
        context.save()
        context.set_source_rgb(.9, .7, .7)
        # draw the y axis
        dx = max(va) - min(va)
        tx = -min(va)/dx
        xzero = (tx * (1 - 2*border_fraction) + border_fraction) * width
        context.move_to(xzero, 0)
        context.line_to(xzero, height)
        context.stroke()
        # draw the x axis
        dy = max(vb) - min(vb)
        ty = -min(vb)/dy
        yzero = (ty * (1 - 2*border_fraction) + border_fraction) * height
        context.move_to(0, yzero)
        context.line_to(width, yzero)
        context.stroke()
        # stop drawing
        context.restore()
    # draw the connections if requested
    if draw_connections:
        # begin drawing
        context.save()
        context.set_source_rgb(.8, .8, .8)
        for i in range(n):
            for j in range(n):
                if i < j and incidence_matrix[i][j] > 0:
                    x, y = rescaled_a[i], rescaled_b[i]
                    nx = (x * (1 - 2*border_fraction) + border_fraction) * width
                    ny = (y * (1 - 2*border_fraction) + border_fraction) * height
                    context.move_to(nx, ny)
                    x, y = rescaled_a[j], rescaled_b[j]
                    nx = (x * (1 - 2*border_fraction) + border_fraction) * width
                    ny = (y * (1 - 2*border_fraction) + border_fraction) * height
                    context.line_to(nx, ny)
                    context.stroke()
        # stop drawing
        context.restore()
    # draw blue internal vertex dots and then green leaf vertex dots
    dot_radius = 2.0
    context.save()
    leaf_rgba = (0.5, 1.0, 0.5, 0.5)
    internal_rgba = (0.5, 0.5, 1.0, 0.5)
    context.set_source_rgba(*internal_rgba)
    for i, (x, y) in enumerate(zip(rescaled_a, rescaled_b)[nleaves:]):
        nx = (x * (1 - 2*border_fraction) + border_fraction) * width
        ny = (y * (1 - 2*border_fraction) + border_fraction) * height
        context.move_to(nx, ny)
        context.arc(nx, ny, dot_radius, 0, 2*math.pi)
        context.fill()
    context.set_source_rgba(*leaf_rgba)
    for i, (x, y) in enumerate(zip(rescaled_a, rescaled_b)[:nleaves]):
        nx = (x * (1 - 2*border_fraction) + border_fraction) * width
        ny = (y * (1 - 2*border_fraction) + border_fraction) * height
        context.move_to(nx, ny)
        context.arc(nx, ny, dot_radius, 0, 2*math.pi)
        context.fill()
    context.restore()
    # draw a scatter plot of the states using the eigenvectors as axes
    for i, (x, y) in enumerate(zip(rescaled_a, rescaled_b)):
        state_string = ordered_names[i]
        nx = (x * (1 - 2*border_fraction) + border_fraction) * width
        ny = (y * (1 - 2*border_fraction) + border_fraction) * height
        context.move_to(nx, ny)
        context.show_text(state_string)
    # get the image string
    return cairo_helper.get_image_string()

def do_projection(D_full, nleaves):
    """
    Compute all projected points onto the plane defined by MDS of the leaves.
    MDS is multidimensional scaling,
    and in this case we are interested in two dimensions because that's what shows best on the monitor.
    @param D_full: the distance matrix as a numpy array relating all vertices including internal vertices
    @param nleaves: the first few indices in D_full represent leaves
    @return: a numpy array where each row is a vertex of the tree viewed as a 2d point
    """
    # Get the points such that the n rows in X are points in n-1 dimensional space.
    X = Euclid.edm_to_points(D_full)
    # Translate all of the points so that the origin is at the centroid of the leaves.
    X -= np.mean(X[:nleaves], 0)
    # Extract the subset of points that define the leaves.
    L = X[:nleaves]
    # Find the orthogonal transformation of the leaves onto their MDS axes.
    # According to the python svd documentation, singular values are sorted most important to least important.
    U, s, Vt = np.linalg.svd(L)
    # Transform all of the points (including the internal vertices) according to this orthogonal transformation.
    # The axes are now the principal axes of the Steiner circumscribed ellipsoid of the leaf vertices.
    # I am using M.T[:2].T to get the first two columns of M.
    vertices_on_plane = np.dot(X, Vt.T).T[:2].T
    return vertices_on_plane

def do_full_projection(D_full, nleaves):
    """
    Compute all projected points onto the subspace defined by the leaves.
    @param D_full: the distance matrix as a numpy array relating all vertices including internal vertices
    @param nleaves: the first few indices in D_full represent leaves
    @return: a numpy array where each row is a vertex of the tree viewed as a point
    """
    # Get the points such that the n rows in X are points in n-1 dimensional space.
    X = Euclid.edm_to_points(D_full)
    # Translate all of the points so that the origin is at the centroid of the leaves.
    X -= np.mean(X[:nleaves], 0)
    # Extract the subset of points that define the leaves.
    L = X[:nleaves]
    # Find the orthogonal transformation of the leaves onto their MDS axes.
    # According to the python svd documentation, singular values are sorted most important to least important.
    U, s, Vt = np.linalg.svd(L)
    # Transform all of the points (including the internal vertices) according to this orthogonal transformation.
    # The axes are now the principal axes of the Steiner circumscribed ellipsoid of the leaf vertices.
    # I am using M.T[:k].T to get the first k columns of M.
    vertices_on_plane = np.dot(X, Vt.T).T[:(nleaves-1)].T
    return vertices_on_plane

def get_response_content(fs):
    # start writing the response type
    response_headers = []
    # read the tree
    tree = NewickIO.parse(fs.tree, FelTree.NewickTree)
    # get the ordered ids and ordered names of the nodes in the tree, starting with the tips
    ordered_name_id_pairs = []
    for node in tree.preorder():
        if node.is_tip():
            name = node.get_name()
            ordered_name_id_pairs.append((name, id(node)))
    for node in tree.preorder():
        if not node.is_tip():
            ordered_name_id_pairs.append(('', id(node)))
    ordered_ids = [id_ for name, id_ in ordered_name_id_pairs]
    ordered_names = [name for name, id_ in ordered_name_id_pairs]
    id_to_index = dict((id_, i) for i, id_ in enumerate(ordered_ids))
    # get the incidence matrix for drawing lines
    n = len(ordered_ids)
    incidence_matrix = [[0]*n for i in range(n)]
    for node in tree.preorder():
        for child in node.gen_children():
            parent_id = id_to_index[id(node)]
            child_id = id_to_index[id(child)]
            incidence_matrix[parent_id][child_id] = 1
            incidence_matrix[child_id][parent_id] = 1
    # get the full distance matrix with ordered indices
    D_full = np.array(tree.get_full_distance_matrix(ordered_ids))
    # get the number of leaves in the tree
    nleaves = len(list(tree.gen_tips()))
    # Compute the projection of all points
    # onto the 2D plane defined by MDS of the leaves.
    projected_points = do_projection(D_full, nleaves)
    # draw the image
    try:
        ext = Form.g_imageformat_to_ext[fs.imageformat]
        image_size = (640, 480)
        return get_image(projected_points, nleaves,
                incidence_matrix, ordered_names, image_size, ext)
    except CairoUtil.CairoUtilError as e:
        raise HandlingError(e)

def examine_projected_distance_matrix():
    tree = NewickIO.parse(g_tree_string, FelTree.NewickTree)
    # get the ordered ids and ordered names of the nodes in the tree, starting with the tips
    ordered_name_id_pairs = []
    for node in tree.preorder():
        if node.is_tip():
            name = node.get_name()
            ordered_name_id_pairs.append((name, id(node)))
    for node in tree.preorder():
        if not node.is_tip():
            ordered_name_id_pairs.append(('', id(node)))
    ordered_ids = [id_ for name, id_ in ordered_name_id_pairs]
    ordered_names = [name for name, id_ in ordered_name_id_pairs]
    id_to_index = dict((id_, i) for i, id_ in enumerate(ordered_ids))
    # get the incidence matrix for drawing lines
    n = len(ordered_ids)
    incidence_matrix = [[0]*n for i in range(n)]
    for node in tree.preorder():
        for child in node.gen_children():
            parent_id = id_to_index[id(node)]
            child_id = id_to_index[id(child)]
            incidence_matrix[parent_id][child_id] = 1
            incidence_matrix[child_id][parent_id] = 1
    # get the full distance matrix with ordered indices
    D_full = np.array(tree.get_full_distance_matrix(ordered_ids))
    # get the number of leaves in the tree
    nleaves = len(list(tree.gen_tips()))
    # compute the projection of all points onto the subspace spanned by the leaves
    projected_points = do_full_projection(D_full, nleaves)
    # compute squared Euclidean distances among the projected points
    D_reconstructed = np.array([[np.dot(b-a, b-a) for b in projected_points] for a in projected_points])
    # set numpy print options for my wide monitor
    np.set_printoptions(linewidth=200)
    # show the original full distance matrix
    print 'original full distance matrix:'
    print np.array(D_full)
    print
    # show the incidence matrix
    print 'incidence matrix:'
    print np.array(incidence_matrix)
    print
    # show the distance matrix constructed from the projected points
    print 'distance matrix reconstructed from projected points:'
    print np.array(D_reconstructed)
    print

def get_ordered_ids(tree):
    """
    Maybe I could use postorder here instead.
    @param tree: a tree
    @return: a list of ids beginning with the leaves
    """
    ordered_ids = []
    ordered_ids.extend(id(node) for node in tree.gen_tips())
    ordered_ids.extend(id(node) for node in tree.gen_internal_nodes())
    return ordered_ids

def examine_mds_splits():
    """
    Examine properties of the hyperplane orthogonal to the MDS axis of a hyperellipse.
    The hyperellipse is the Steiner circumscribed hyperellipse that intersects
    points of the embedded leaves of a tree.
    Earlier results show that the hyperplane orthogonal to the principal
    axis of this hyperellipse should separate the leaves in a way that is compatible
    with the topology of the tree.
    Here we investigate the conjecture that this same hyperplane
    also splits internal vertices in a way that is compatible with the topology of the tree.
    """
    count = 0
    ncontrol_noneuclidean_counterexamples = 0
    ncontrol_secondary_counterexamples = 0
    print 'Does the principal hyperplane of the leaves always intersect the tree at exactly one point?'
    print 'Press control-C to stop looking for a counterexample...'
    try:
        while True:
            # pick a random number of taxa to use as leaves in the tree
            ntaxa = random.randrange(3, 12)
            # sample an xtree with exponentially distributed branch lengths
            xtree = TreeSampler.sample_agglomerated_tree(ntaxa)
            for branch in xtree.get_branches():
                mu = 2.0
                branch.length = random.expovariate(1/mu)
            # convert the xtree to a FelTree so we can use the internal vertices
            tree_string = xtree.get_newick_string()
            tree = NewickIO.parse(tree_string, FelTree.NewickTree)
            # get the full id splits of the tree, including internal nodes
            id_set = set(id(node) for node in tree.preorder())
            d = TreeComparison._get_branch_id_to_node_id_set(tree)
            full_id_splits = set(frozenset((frozenset(x), frozenset(id_set-x))) for x in d.values())
            # get ordered ids and the number of leaves
            ordered_ids = get_ordered_ids(tree)
            nleaves = len(list(tree.gen_tips()))
            # get the projection
            D_full = np.array(tree.get_full_distance_matrix(ordered_ids))
            projected_points = do_projection(D_full, nleaves)
            # get the split implied by the principal hyperplane of the leaves
            left_ids = set(i for i, point in zip(ordered_ids, projected_points) if point[0] < 0)
            right_ids = id_set - left_ids
            split = frozenset((frozenset(left_ids), frozenset(right_ids)))
            # if the split is not compatible with the tree then we have found a counterexample
            if split not in full_id_splits:
                print 'counterexample:'
                print tree_string
                break
            # now do a control where I look at the wrong eigenvector
            left_ids = set(i for i, point in zip(ordered_ids, projected_points) if point[1] < 0)
            right_ids = id_set - left_ids
            split = frozenset((frozenset(left_ids), frozenset(right_ids)))
            if split not in full_id_splits:
                ncontrol_secondary_counterexamples += 1
            # now do a control that should provide the occasional counterexample
            D_control = np.sqrt(D_full)
            projected_points = do_projection(D_control, nleaves)
            left_ids = set(i for i, point in zip(ordered_ids, projected_points) if point[0] < 0)
            right_ids = id_set - left_ids
            split = frozenset((frozenset(left_ids), frozenset(right_ids)))
            if split not in full_id_splits:
                ncontrol_noneuclidean_counterexamples += 1
            # increment the count
            count += 1
    except KeyboardInterrupt, e:
        print 'Checked', count, 'trees and found no counterexample.'
        print 'Found', ncontrol_secondary_counterexamples, 'control counterexamples where I use the wrong eigenvector.'
        print 'Found', ncontrol_noneuclidean_counterexamples, 'control counterexamples where I use the wrong distance matrix.'

def main():
    #examine_projected_distance_matrix()
    examine_mds_splits()

if __name__ == '__main__':
    main()
