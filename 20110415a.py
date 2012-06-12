"""
Animate a 2D MDS tree as a branch length changes.

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
import argparse

import numpy as np
import cairo

import Form
import FormOut
import Ftree
import FtreeIO
import Euclid
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
            Form.Float('scale', 'scale the image of the tree by this factor',
                200.0, low_exclusive=0.0),
            Form.Float('frame_progress',
                'animation frame progress between 0.0 and 1.0',
                0.5, low_inclusive=0.0, high_inclusive=1.0),
            Form.SingleLine('branch_name',
                'name of the vertex associated with the variable branch', '1'),
            Form.Float('final_length',
                'final branch length', 10.0, low_exclusive=0.0),
            Form.Integer('x_axis',
                'x axis projection (1 is Fiedler)', 1, low=1),
            Form.Integer('y_axis',
                'y axis projection (1 is Fiedler)', 2, low=1),
            Form.ImageFormat()]
    return form_objects

def get_form_out():
    return FormOut.Image('frame')

def reflect_to_match(A, B):
    """
    The idea is to make a new matrix by multiplying columns of A by -1 or 1.
    In this new matrix each column vector will be pointing in the
    same approximate direction as the corresponding column
    vector in B.
    @param A: the new matrix
    @param B: the old matrix
    @return: a reflection of the new matrix
    """
    return np.dot(A, np.diag(np.sign(np.diag(np.dot(A.T, B)))))

def raw_sigmoid(x, alpha):
    """
    @param x: something between -0.5 and 0.5
    @param alpha: larger alpha means a steeper shift
    @return: something small but monotonic
    """
    return 1.0 / (1.0 + math.exp(-alpha*x)) - 0.5

def sigmoid(x, alpha=12):
    """
    @param x: something between 0 and 1
    @param alpha: larger alpha means a steeper shift
    @return: something between 0 and 1
    """
    height = raw_sigmoid(0.5, alpha) - raw_sigmoid(-0.5, alpha)
    t = raw_sigmoid(x - 0.5, alpha) / height + 0.5
    if t < 0 or t > 1:
        raise ValueError('oops i screwed up the response curve')
    return t

def get_edge(R, N, name):
    pairs = [(a,b) for a, b in R if N.get(b, None) == name]
    if len(pairs) > 1:
        raise ValueError(
                'expected the vertex to define a single variable branch')
    if len(pairs) < 1:
        raise ValueError('the provided vertex is not associated with a branch')
    return frozenset(pairs[0])

def get_response_content(fs):
    # define the requested physical size of the images (in pixels)
    physical_size = (640, 480)
    # get the directed edges and the branch lengths and vertex names
    R, B, N = FtreeIO.newick_to_RBN(fs.tree_string)
    # get the requested undirected edge
    edge = get_edge(R, N, fs.branch_name)
    # get the undirected tree topology
    T = Ftree.R_to_T(R)
    # get the leaves and the vertices of articulation
    leaves = Ftree.T_to_leaves(T)
    internal = Ftree.T_to_internal_vertices(T)
    vertices = leaves + internal
    nleaves = len(leaves)
    v_to_index = Ftree.invseq(vertices)
    # get the requested indices
    x_index = fs.x_axis - 1
    y_index = fs.y_axis - 1
    if x_index >= nleaves-1 or y_index >= nleaves-1:
        raise ValueError(
                'projection indices must be smaller than the number of leaves')
    # adjust the branch length
    initial_length = B[edge]
    t = sigmoid(fs.frame_progress)
    B[edge] = (1-t)*initial_length + t*fs.final_length
    # get the points
    w, v = Ftree.TB_to_harmonic_extension(T, B, leaves, internal)
    X_full = np.dot(v, np.diag(np.reciprocal(np.sqrt(w))))
    X = np.vstack([X_full[:,x_index], X_full[:,y_index]]).T
    # draw the image
    ext = Form.g_imageformat_to_ext[fs.imageformat]
    return get_animation_frame(ext, physical_size, fs.scale, v_to_index, T, X)

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

def main(args):
    # do some validation
    if args.nframes < 2:
        raise ValueError('nframes should be at least 2')
    # define the requested physical size of the images (in pixels)
    physical_size = (args.physical_width, args.physical_height)
    # get the directed edges and the branch lengths and vertex names
    R, B, N = FtreeIO.newick_to_RBN(args.tree)
    # get the requested undirected edge
    edge = get_edge(R, N, args.branch_name)
    initial_length = B[edge]
    # get the undirected tree topology
    T = Ftree.R_to_T(R)
    # get the leaves and the vertices of articulation
    leaves = Ftree.T_to_leaves(T)
    internal = Ftree.T_to_internal_vertices(T)
    vertices = leaves + internal
    nleaves = len(leaves)
    v_to_index = Ftree.invseq(vertices)
    # get the requested indices
    x_index = args.x_axis - 1
    y_index = args.y_axis - 1
    if x_index >= nleaves-1 or y_index >= nleaves-1:
        raise ValueError(
                'projection indices must be smaller than the number of leaves')
    X_prev = None
    # create the animation frames and write them as image files
    pbar = Progress.Bar(args.nframes)
    for frame_index in range(args.nframes):
        linear_progress = frame_index / float(args.nframes - 1)
        if args.interpolation == 'sigmoid':
            t = sigmoid(linear_progress)
        else:
            t = linear_progress
        B[edge] = (1-t)*initial_length + t*args.final_length
        w, v = Ftree.TB_to_harmonic_extension(T, B, leaves, internal)
        X_full = np.dot(v, np.diag(np.reciprocal(np.sqrt(w))))
        X = np.vstack([X_full[:,x_index], X_full[:,y_index]]).T
        if X_prev is not None:
            X = reflect_to_match(X, X_prev)
        X_prev = X
        image_string = get_animation_frame(
                args.image_format, physical_size, args.scale,
                v_to_index, T, X)
        image_filename = 'frame-%04d.%s' % (frame_index, args.image_format)
        image_pathname = os.path.join(args.output_directory, image_filename)
        with open(image_pathname, 'wb') as fout:
            fout.write(image_string)
        pbar.update(frame_index+1)
    pbar.finish()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scale', type=float, default=200.0,
            help='define the drawing scale') 
    parser.add_argument('--physical_width', type=int, default=480,
            help='width (pixels)') 
    parser.add_argument('--physical_height', type=int, default=360,
            help='height (pixels)') 
    parser.add_argument('--tree', default=g_tree_string,
            help='newick tree with branch lengths')
    parser.add_argument('--branch_name', default='1',
            help='name of the vertex associated with the variable branch')
    parser.add_argument('--final_length', type=float, default=10.0,
            help='final branch length')
    parser.add_argument('--x_axis', type=int, default=1,
            help='x axis projection (1 is Fiedler)')
    parser.add_argument('--y_axis', type=int, default=2,
            help='y axis projection (1 is Fiedler)')
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

