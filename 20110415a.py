"""Animate a 2D MDS tree as a branch length changes. [UNFINISHED]

Create a tree MDS animation
showing a progressive branch length change.
A sequence of .png files should be written
to some existing specified output directory.
If a web interface is used, maybe show one frame
at some stage between 0 and 1.
Input for web usage:
the tree, a progress fraction, a scaling factor,
image format, and image delivery.
Input for command line usage:
the path to the output directory for the images, a scaling factor, a tree,
and a physical width and height.
To convert a sequence of png images to an mpeg video:
ffmpeg -i frames/frame-%04d.png test.mpg
Resolutions preferred by YouTube are 640x360 and 480x360.
"""


from StringIO import StringIO
import os
import math

import numpy as np
import scipy
import cairo
import argparse

from SnippetUtil import HandlingError
import SnippetUtil
import Form
import FormOut
import Ftree
import FtreeIO
import Euclid
import CairoUtil
import Progress
import MatrixUtil
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
            Form.Float('scale', 'scale the image of the tree by this factor',
                200.0, low_exclusive=0.0),
            Form.Float('progress', 'animation progress between 0.0 and 1.0',
                0.5, low_inclusive=0.0, high_inclusive=1.0),
            Form.SingleLine('branch_name',
                'name of the vertex associated with the variable branch', '1'),
            Form.Float('final_blen',
                'final branch length', 10.0, low_exclusive=0.0),
            Form.Integer('x_axis',
                'x axis projection (1 is Fiedler)', 1, low=1),
            Form.Integer('y_axis',
                'y axis projection (1 is Fiedler)', 2, low=1),
            Form.ImageFormat(),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.Image('frame')

def get_response_content(fs):
    # define the requested physical size of the images (in pixels)
    physical_size = (640, 480)
    # get the directed edges and the branch lengths and vertex names
    R, B, N = FtreeIO.newick_to_RBN(fs.tree_string)
    # get the requested undirected edge
    pairs = [(a,b) for a, b in R if N.get(b, None) == fs.branch_name]
    if len(pairs) > 1:
        msg = 'expected the vertex to define a single variable branch'
        raise ValueError(msg)
    if len(pairs) < 1:
        msg = 'the provided vertex is not associated with a branch'
        raise ValueError(msg)
    edge = frozenset(pairs[0])
    # get the undirected tree topology
    T = Ftree.R_to_T(R)
    # get the leaves and the vertices of articulation
    leaves = Ftree.T_to_leaves(T)
    internal = Ftree.T_to_internal_vertices(T)
    vertices = leaves + internal
    nleaves = len(leaves)
    nvertices = len(vertices)
    v_to_index = Ftree.invseq(vertices)
    # get the requested indices
    x_index = fs.x_axis - 1
    y_index = fs.y_axis - 1
    if x_index >= nleaves-1 or y_index >= nleaves-1:
        msg = 'projection indices must be smaller than the number of leaves'
        raise ValueError(msg)
    w, v = Ftree.TB_to_harmonic_extension(T, B, leaves, internal)
    print T
    print B
    print w
    print v
    X_full = np.dot(v, np.diag(np.reciprocal(np.sqrt(w))))
    X = np.vstack([X_full[:,x_index], X_full[:,y_index]]).T
    # draw the image
    ext = Form.g_imageformat_to_ext[fs.imageformat]
    return get_animation_frame(ext, physical_size, fs.scale,
            v_to_index, T, X)

def get_animation_frame(
        image_format, physical_size, scale,
        v_to_index, T, X):
    """
    This function is about drawing the tree.
    @param image_format: the image extension
    @param physical_size: the width and height of the image in pixels
    @param scale: a scaling factor
    @param v_to_index: maps vertices to their index
    @param T: defines the connectivity of the tree
    @param X: an array of 2D points
    @return: the animation frame as an image as a string
    """
    # before we begin drawing we need to create the cairo surface and context
    cairo_helper = CairoUtil.CairoHelper(image_format)
    surface = cairo_helper.create_surface(physical_size[0], physical_size[1])
    context = cairo.Context(surface)
    # define some helper variables
    x0 = physical_size[0] / 2.0
    y0 = physical_size[1] / 2.0
    npoints = len(X)
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
    for va, vb in T:
        a = v_to_index[va]
        b = v_to_index[vb]
        ax, ay = X[a].tolist()
        bx, by = X[b].tolist()
        context.move_to(x0 + ax*scale, y0 + ay*scale)
        context.line_to(x0 + bx*scale, y0 + by*scale)
        context.stroke()
    context.restore()
    # create the image
    return cairo_helper.get_image_string()

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
    # Get ordered ids with the leaves first,
    # and get the corresponding distance matrix.
    ordered_ids = get_ordered_ids(tree)
    D = np.array(tree.get_partial_distance_matrix(ordered_ids))
    index_edges = get_index_edges(tree, ordered_ids)
    # Create the reference points
    # so that the video frames are not reflected arbitrarily.
    reference_points = Euclid.edm_to_points(D).T[:3].T
    # create the animation frames and write them as image files
    pbar = Progress.Bar(args.nframes)
    for frame_index in range(args.nframes):
        linear_progress = frame_index / float(args.nframes - 1)
        if args.interpolation == 'sigmoid':
            progress = sigmoid(linear_progress)
        else:
            progress = linear_progress
        mass_vector = get_mass_vector(nvertices, nleaves, progress)
        points = get_canonical_3d_mds(D, mass_vector, reference_points)
        crossings = get_crossings(index_edges, points)
        image_string = get_animation_frame(
                args.image_format, physical_size, args.scale,
                mass_vector, index_edges, points, crossings)
        image_filename = 'frame-%04d.%s' % (frame_index, args.image_format)
        image_pathname = os.path.join(args.output_directory, image_filename)
        with open(image_pathname, 'wb') as fout:
            fout.write(image_string)
        pbar.update(frame_index+1)
    pbar.finish()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scale', type=float, default=1.0,
            help='define the drawing scale') 
    parser.add_argument('--physical_width', type=int, default=480,
            help='width (pixels)') 
    parser.add_argument('--physical_height', type=int, default=360,
            help='height (pixels)') 
    parser.add_argument('--tree', default=g_tree_string,
            help='newick tree with branch lengths')
    parser.add_argument('--image_format', default='png',
            choices=('png', 'svg', 'ps', 'pdf'),
            help='image format')
    parser.add_argument('--nframes', type=int, default=100,
            help='number of animation frames (image files) to create') 
    parser.add_argument('--interpolation', default='sigmoid',
            choices=('sigmoid', 'linear'),
            help='weights change according to this function')
    parser.add_argument('output_directory',
            help='path to the output directory for .png frames')
    main(parser.parse_args())

