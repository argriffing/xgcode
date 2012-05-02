"""
Animate a 3D tree.

To convert a sequence of png images to an mpeg video:
ffmpeg -i frames/frame-%04d.png test.mpg
Resolutions preferred by YouTube are 640x360 and 480x360. 
"""


from StringIO import StringIO
import os
import math
import argparse

import numpy as np
import scipy.stats
import cairo

import Form
import FormOut
import Ftree
import FtreeIO
import Euclid
import CairoUtil
import Progress
import const
import TreeProjection
import NewickIO

#g_tree_string = const.read('20100730g').rstrip()
g_tree_string = NewickIO.daylight_example_tree

g_default_yaw = 0
g_default_pitch = 0.2
g_default_eigenvector_index = 2

def get_form():
    """
    @return: a list of form objects
    """
    # define the form objects
    form_objects = [
            Form.MultiLine('tree_string', 'newick tree',
                g_tree_string),
            #Form.Float('scale', 'scale the image of the tree by this factor',
            #200.0, low_exclusive=0.0),
            Form.Integer('eigenvector_index',
                'eigenvector index (1 is Fiedler)',
                g_default_eigenvector_index, low=1),
            Form.Float('yaw', 'yaw', g_default_yaw),
            Form.Float('pitch', 'pitch', g_default_pitch),
            Form.ImageFormat()]
    return form_objects

def get_form_out():
    return FormOut.Image('frame')

def get_response_content(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: the response
    """
    physical_size = (640, 480)
    scale = 10.0
    return get_animation_frame(
            fs.imageformat, physical_size, scale,
            fs.tree_string, fs.eigenvector_index, fs.yaw, fs.pitch)

def t_to_yaw(t):
    """
    @param t: animation progress between 0 and 1
    @return: yaw
    """
    yaw = math.pi*2*t
    return yaw

def t_to_pitch(t):
    """
    @param t: animation progress between 0 and 1
    @return: pitch
    """
    pitch = (math.pi / 4.0) * math.sin(math.pi * 2.0 * 3 * t)
    return pitch

def get_animation_frame(
        image_format, physical_size, scale,
        newick, eigenvector_index, yaw, pitch):
    """
    This function is about drawing the tree.
    @param image_format: the image extension
    @param physical_size: the width and height of the image in pixels
    @param scale: a scaling factor
    @return: the animation frame as an image as a string
    """
    # before we begin drawing we need to create the cairo surface and context
    cairo_helper = CairoUtil.CairoHelper(image_format)
    surface = cairo_helper.create_surface(physical_size[0], physical_size[1])
    ctx = cairo.Context(surface)
    # draw a white background 
    ctx.save() 
    ctx.set_source_rgb(1, 1, 1)
    ctx.paint() 
    ctx.restore() 
    # define some helper variables
    x0 = physical_size[0] / 2.0
    y0 = physical_size[1] / 2.0
    # translate to the center of the frame
    ctx.translate(x0, y0)
    ctx.scale(1, -1)
    # draw the info
    TreeProjection.draw_cairo_frame(
            ctx, scale, newick, eigenvector_index, yaw, pitch)
    # create the image
    return cairo_helper.get_image_string()

def main(args):
    # do some validation
    if args.nframes < 2:
        raise ValueError('nframes should be at least 2')
    # define the requested physical size of the images (in pixels)
    physical_size = (args.physical_width, args.physical_height)
    # create the animation frames and write them as image files
    pbar = Progress.Bar(args.nframes)
    for frame_index in range(args.nframes):
        t = frame_index / float(args.nframes - 1)
        image_string = get_animation_frame(
                args.image_format, physical_size, args.scale,
                args.tree, args.eigenvector_index,
                t_to_yaw(t), t_to_pitch(t))
        image_filename = 'frame-%04d.%s' % (frame_index, args.image_format)
        image_pathname = os.path.join(args.output_directory, image_filename)
        with open(image_pathname, 'wb') as fout:
            fout.write(image_string)
        pbar.update(frame_index+1)
    pbar.finish()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scale', type=float, default=10.0,
            help='define the drawing scale') 
    parser.add_argument('--physical_width', type=int, default=480,
            help='width (pixels)') 
    parser.add_argument('--physical_height', type=int, default=360,
            help='height (pixels)') 
    parser.add_argument('--tree', default=g_tree_string,
            help='newick tree with branch lengths')
    parser.add_argument('--eigenvector_index',
            default=g_default_eigenvector_index, type=int,
            help='eigenvector index (1 is Fiedler)')
    parser.add_argument('--image_format', default='png',
            choices=('png', 'svg', 'ps', 'pdf'),
            help='image format')
    parser.add_argument('--nframes', type=int, default=100,
            help='number of animation frames (image files) to create') 
    parser.add_argument('output_directory',
            help='path to the output directory for .png frames')
    main(parser.parse_args())

