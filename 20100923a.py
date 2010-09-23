"""Plot an MDS spectrogram showing eigenvalues as tree tips are upweighted. [UNFINISHED]
"""


from StringIO import StringIO
import random
import os
import math
from itertools import product

import numpy as np
import cairo
import argparse

from SnippetUtil import HandlingError
import SnippetUtil
import Form
import FormOut
import NewickIO
import Euclid
import NeighborJoining
import FelTree
import CairoUtil
import Progress
import const

g_tree_string = const.read('20100730g').rstrip()


def get_form():
    """
    @return: a list of form objects
    """
    # define the form objects
    form_objects = [
            Form.MultiLine('tree_string', 'newick tree',
                g_tree_string),
            Form.ImageFormat(),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.Image('spectrogram')

def create_test_image(image_format, width, height):
    """
    @param image_format: a string like an image extension
    @param width: width in pixels
    @param height: height in pixels
    @return: the contents of the image as a string
    """
    # do some initialization
    cairo_helper = CairoUtil.CairoHelper(image_format)
    surface = cairo_helper.create_surface(width, height)
    context = cairo.Context(surface)
    # draw an off-white background
    context.save()
    context.set_source_rgb(.9, .9, .9)
    context.paint()
    context.restore()
    # get the contents of the image
    return cairo_helper.get_image_string()

def get_response_content(fs):
    # define the requested physical size of the images (in pixels)
    physical_size = (640, 480)
    # build the newick tree from the string
    tree = NewickIO.parse(fs.tree_string, FelTree.NewickTree)
    nvertices = len(list(tree.preorder()))
    nleaves = len(list(tree.gen_tips()))
    # Get ordered ids with the leaves first,
    # and get the corresponding distance matrix.
    ordered_ids = get_ordered_ids(tree)
    D = np.array(tree.get_partial_distance_matrix(ordered_ids))
    index_edges = get_index_edges(tree, ordered_ids)
    # Create the reference points so that the video frames
    # are not reflected arbitrarily.
    reference_points = Euclid.edm_to_points(D).T[:2].T
    # draw the image
    ext = Form.g_imageformat_to_ext[fs.imageformat]
    mass_vector = get_mass_vector(nvertices, nleaves, fs.progress)
    points = get_canonical_2d_mds(D, mass_vector, reference_points)
    return get_animation_frame(ext, physical_size, fs.scale,
            mass_vector, index_edges, points)

def reflect_to_reference(points, reference_points):
    """
    For 2D points, try each of the 4 possible combinations reflections across the axes.
    Use the reflection that gives points closest to the reference points.
    @param points: rows are 2D points
    @param reference_points: rows are 2D points
    @return: 
    """
    if points.shape != reference_points.shape:
        raise ValueError('the point array and the reference point array should have the same shape')
    if len(points.shape) != 2:
        raise ValueError('the points argument should be a matrix-like numpy array')
    if points.shape[1] != 2:
        raise ValueError('the points should be in 2D space')
    reflectors = np.array(list(product((-1,1), repeat=2)))
    best_error, best_reflector = min((np.linalg.norm(points*r - reference_points), r) for r in reflectors)
    return points * best_reflector

def get_index_edges(tree, ordered_ids):
    """
    Given a tree and some ordered ids, get a collection of edges defined on indices.
    @param tree: the tree object
    @param ordered_ids: the returned index pairs are with respect to this sequence
    @return: a collection of index pairs defining edges
    """
    # map ids to indices
    id_to_index = dict((myid, index) for index, myid in enumerate(ordered_ids))
    # each edge in this set is a frozenset of two indices
    index_edges = set()
    for node in tree.preorder():
        index = id_to_index[id(node)]
        for neighbor in node.gen_neighbors():
            neighbor_index = id_to_index[id(neighbor)] 
            index_edges.add(frozenset([index, neighbor_index]))
    return index_edges

def get_mass_vector(nvertices, nleaves, progress):
    """
    @param nvertices: the number of vertices in the tree, including internal vertices
    @param nleaves: the number of leaves in the tree
    @param progress: between 0.0 (uniform weights) and 1.0 (uniform tip weights, zero internal weights)
    @return: a nonnegative mass vector that sums to 1.0 for weighting the MDS
    """
    # tips get weights proportional to 1.0 and internal vertices get weights proportional to (1-progress)
    mass_vector = np.ones(nvertices, dtype=float)
    for i in range(nleaves, nvertices):
        mass_vector[i] = 1-progress
    return mass_vector / sum(mass_vector)

def get_canonical_2d_mds(D, m, reference_points):
    """
    This function is about projecting the points.
    It is like MDS except the reflections across the axes are not arbitrary.
    Also it only uses the first two axes.
    @param D: the full distance matrix
    @param m: the mass vector
    @param reference_points: a 2D reference projection of vertices of the tree
    @return: the weighted MDS points as a numpy matrix
    """
    X = Euclid.edm_to_weighted_points(D, m)
    return reflect_to_reference(X.T[:2].T, reference_points)

def get_animation_frame(image_format, physical_size, scale, mass_vector, index_edges, points):
    """
    This function is about drawing the tree.
    @param image_format: the image extension
    @param physical_size: the width and height of the image in pixels
    @param scale: a scaling factor
    @param mass_vector: use this for visualizing the weights of the vertices
    @param index_edges: defines the connectivity of the tree
    @param points: an array of 2D points, the first few of which are leaves
    @return: the animation frame as an image as a string
    """
    # before we begin drawing we need to create the cairo surface and context
    cairo_helper = CairoUtil.CairoHelper(image_format)
    surface = cairo_helper.create_surface(physical_size[0], physical_size[1])
    context = cairo.Context(surface)
    # define some helper variables
    x0 = physical_size[0] / 2.0
    y0 = physical_size[1] / 2.0
    npoints = len(points)
    # draw an off-white background
    context.save()
    context.set_source_rgb(.9, .9, .9)
    context.paint()
    context.restore()
    # draw the axes which are always in the center of the image
    context.save()
    context.set_source_rgb(.9, .7, .7)
    context.move_to(x0, 0)
    context.line_to(x0, physical_size[1])
    context.stroke()
    context.move_to(0, y0)
    context.line_to(physical_size[0], y0)
    context.stroke()
    context.restore()
    # draw the edges
    context.save()
    context.set_source_rgb(.8, .8, .8)
    for edge in index_edges:
        ai, bi = tuple(edge)
        ax, ay = points[ai].tolist()
        bx, by = points[bi].tolist()
        context.move_to(x0 + ax*scale, y0 + ay*scale)
        context.line_to(x0 + bx*scale, y0 + by*scale)
        context.stroke()
    context.restore()
    # draw vertices as translucent circles with radius defined by the mass vector
    context.save()
    context.set_source_rgba(0.2, 0.2, 1.0, 0.5)
    for point, mass in zip(points, mass_vector):
        if mass:
            x, y = point.tolist()
            nx = x0 + x*scale
            ny = y0 + y*scale
            dot_radius = 3*mass*npoints
            context.arc(nx, ny, dot_radius, 0, 2*math.pi)
            context.fill()
    context.restore()
    # create the image
    return cairo_helper.get_image_string()

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

def sigmoid(x):
    t = (x - .5) * 12
    return 1.0 / (1.0 + math.exp(-t))

def main(args):
    # do some validation
    if args.nframes < 2:
        raise ValueError('nframes should be at least 2')
    # define the requested physical size of the images (in pixels)
    physical_size = (args.physical_width, args.physical_height)
    # build the newick tree from the string
    tree = NewickIO.parse(args.tree, FelTree.NewickTree)
    nvertices = len(list(tree.preorder()))
    nleaves = len(list(tree.gen_tips()))
    # get ordered ids with the leaves first, and get the corresponding distance matrix
    ordered_ids = get_ordered_ids(tree)
    D = np.array(tree.get_partial_distance_matrix(ordered_ids))
    index_edges = get_index_edges(tree, ordered_ids)
    # create the reference points so that the video frames are not reflected arbitrarily
    reference_points = Euclid.edm_to_points(D).T[:2].T
    # create the animation frames and write them as image files
    pbar = Progress.Bar(args.nframes)
    for frame_index in range(args.nframes):
        linear_progress = frame_index / float(args.nframes - 1)
        if args.interpolation == 'sigmoid':
            progress = sigmoid(linear_progress)
        else:
            progress = linear_progress
        mass_vector = get_mass_vector(nvertices, nleaves, progress)
        points = get_canonical_2d_mds(D, mass_vector, reference_points)
        image_string = get_animation_frame(args.image_format, physical_size, args.scale, mass_vector, index_edges, points)
        image_filename = 'frame-%04d.%s' % (frame_index, args.image_format)
        image_pathname = os.path.join(args.output_directory, image_filename)
        with open(image_pathname, 'wb') as fout:
            fout.write(image_string)
        pbar.update(frame_index+1)
    pbar.finish()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scale', type=float, default=1.0, help='define the drawing scale') 
    parser.add_argument('--physical_width', type=int, default=480, help='width (pixels)') 
    parser.add_argument('--physical_height', type=int, default=360, help='height (pixels)') 
    parser.add_argument('--tree', default=g_tree_string, help='newick tree with branch lengths')
    parser.add_argument('--image_format', default='png', choices=('png', 'svg', 'ps', 'pdf'), help='image format')
    parser.add_argument('--nframes', type=int, default=100, help='number of animation frames (image files) to create') 
    parser.add_argument('--interpolation', default='sigmoid', choices=('sigmoid', 'linear'), help='weights change according to this function')
    parser.add_argument('output_directory', help='path to the directory to which the .png frames should be written')
    main(parser.parse_args())

