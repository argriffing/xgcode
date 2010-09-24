"""Plot an MDS spectrogram showing eigenvalues as tree tips are upweighted. [UNFINISHED]

Consider three eigenvalue scalings.
First, scale by the max over all eigenvalues
for all duplication counts.
Second, scale by the duplication count.
Third, scale by the max over all eigenvalues
for each duplication count individually.
"""


from StringIO import StringIO
import os
import math

import numpy as np
import cairo
import argparse
import scipy

from SnippetUtil import HandlingError
import SnippetUtil
import Form
import FormOut
import NewickIO
import Euclid
import FelTree
import CairoUtil
import const
import MatrixUtil

g_tree_string = const.read('20100730g').rstrip()


def get_form():
    """
    @return: a list of form objects
    """
    # define the form objects
    form_objects = [
            Form.MultiLine('tree_string', 'newick tree',
                g_tree_string),
            Form.Integer('border', 'number of border pixels', 5),
            Form.Integer('hticks', 'number of ticks on horizontal axis', 200),
            Form.Integer('denom', 'denominator of weight ratio', 4),
            Form.ImageFormat(),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.Image('spectrogram')

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
    # get the image extension
    ext = Form.g_imageformat_to_ext[fs.imageformat]
    # get the scaling factors and offsets
    if fs.hticks < 2:
        msg = 'expected at least two ticks on the horizontal axis'
        raise HandlingError(msg)
    width, height = physical_size
    xoffset = fs.border
    yoffset = fs.border
    yscale = float(height - 2*fs.border)
    xscale = (width - 2*fs.border) / float(fs.hticks - 1)
    # draw the image
    return create_image(ext, physical_size,
            xscale, yscale, xoffset, yoffset,
            D, nleaves, fs.hticks, fs.denom)

def expand_vdup(vdup):
    """
    @param vdup: a vector of row duplication counts
    @return: a vector of indices
    """
    indices = []
    for i, count in enumerate(vdup):
        indices.extend([i]*count)
    return indices

def get_dup_distance_matrix(D_in, vdup):
    """
    @param D_in: a distance matrix
    @param vdup: the number of multiples of each original row
    @return: a large distance matrix with many row and column duplicates
    """
    D = np.zeros((sum(vdup), sum(vdup)))
    for i_sink, i_source in enumerate(expand_vdup(vdup)):
        for j_sink, j_source in enumerate(expand_vdup(vdup)):
            D[i_sink, j_sink] = D_in[i_source, j_source]
    return D

def get_augmented_spectrum(D_in, ntips, ndup_tip, ndup_internal):
    """
    The tips are the first indices of the original distance matrix.
    @param D_in: the original distance matrix
    @param ntips: the number of tips in the tree
    @param ndup_tip: the total number of repeats per tip node
    @param ndup_internal: the total number of repeats per internal node
    @return: eigenvalues of Gower's centered augmented distance matrix
    """
    vdup = [ndup_tip]*ntips + [ndup_internal]*(len(D_in) - ntips)
    D = get_dup_distance_matrix(D_in, vdup)
    G = -0.5 * MatrixUtil.double_centered(D)
    w = scipy.linalg.eigh(G, eigvals_only=True)
    return (w * ndup_internal) / ndup_tip

def create_image(ext, size,
        xscale, yscale, xoffset, yoffset,
        D_in, ntips, nticks_horz, denom):
    """
    @param ext: image extension
    @param size: width and height of the image in pixels
    @param D_in: the original distance matrix
    @param ntips: the number of tips in the tree
    @param nticks_horz: the number of ticks on the horizontal axis
    @return: image file contents
    """
    # get all of the scaled eigenvalues we want to plot
    ndups_list = range(denom, nticks_horz+denom)
    raw_eigenvalue_grid = [
            get_augmented_spectrum(D_in, ntips, k, denom) for k in ndups_list]
    max_eigenvalue = max(max(ws) for ws in raw_eigenvalue_grid)
    scaled_eigenvalue_grid = [
            np.array(w) / max_eigenvalue for w in raw_eigenvalue_grid]
    # before we begin drawing we need to create the cairo surface and context
    cairo_helper = CairoUtil.CairoHelper(ext)
    surface = cairo_helper.create_surface(size[0], size[1])
    context = cairo.Context(surface)
    # draw an off-white background
    context.save()
    context.set_source_rgb(.9, .9, .9)
    context.paint()
    context.restore()
    # draw eigenvalues as translucent disks
    context.save()
    context.set_source_rgba(0.2, 0.2, 1.0, 0.5)
    dot_radius = 2
    for htick, ws in enumerate(scaled_eigenvalue_grid):
        x = xoffset + htick * xscale
        for w in ws:
            y = size[1] - (yoffset + w * yscale)
            context.arc(x, y, dot_radius, 0, 2*math.pi)
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
    tips = [id(node) for node in tree.gen_tips()]
    internal = [id(node) for node in tree.gen_internal_nodes()]
    return tips + internal
