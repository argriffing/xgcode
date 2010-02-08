"""Plot the first few eigenfunctions of the Laplacian of a path.
"""

from StringIO import StringIO
from itertools import product

import numpy as np
import cairo

from SnippetUtil import HandlingError
import Form
import Euclid
import CairoUtil
import iterutils

g_naxes = 4

def get_form():
    """
    @return: a list of form objects
    """
    # define the form objects
    form_objects = [
            Form.Integer('nvertices', 'use this many vertices',
                20, low=g_naxes+1, high=100),
            Form.RadioGroup('eigenvalue_option', 'eigenvalue usage', [
                Form.RadioItem('eigenvalue_scaling',
                    'scale by square roots of eigenvalues', True),
                Form.RadioItem('eigenvalue_ignoring',
                    'ignore eigenvalues')]),
            Form.RadioGroup('xaxis_option', 'plotting options', [
                Form.RadioItem('xaxis_length',
                    'x axis is the path distance', True),
                Form.RadioItem('xaxis_firstvector',
                    'x axis is the first MDS vector')]),
            Form.RadioGroup('imageformat', 'image format', [
                Form.RadioItem('png', 'png', True),
                Form.RadioItem('svg', 'svg'),
                Form.RadioItem('pdf', 'pdf'),
                Form.RadioItem('ps', 'ps')]),
            Form.RadioGroup('contentdisposition', 'image delivery', [
                Form.RadioItem('inline', 'view the image', True),
                Form.RadioItem('attachment', 'download the image')])]
    return form_objects

def create_laplacian_matrix(nvertices):
    """
    @param affinity: affinity between adjacent vertices
    @param nvertices: the number of vertices in the graph
    @return: a numpy matrix
    """
    affinity = nvertices * 2.0
    A = np.zeros((nvertices, nvertices), dtype=float)
    for i, j in iterutils.pairwise(range(nvertices)):
        A[i,j] = affinity
        A[j,i] = affinity
    L = Euclid.adjacency_to_laplacian(A)
    return L

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # define the requested physical size of the images (in pixels)
    physical_size = (640, 480)
    # get the points
    L = create_laplacian_matrix(fs.nvertices)
    D = Euclid.laplacian_to_edm(L)
    HSH = Euclid.edm_to_dccov(D)
    W, VT = np.linalg.eigh(HSH)
    V = VT.T.tolist()
    if fs.eigenvalue_scaling:
        vectors = [np.array(v)*w for w, v in list(reversed(sorted(zip(np.sqrt(W), V))))[:-1]]
    else:
        vectors = [np.array(v) for w, v in list(reversed(sorted(zip(np.sqrt(W), V))))[:-1]]
    X = np.array(zip(*vectors))
    # transform the points to eigenfunctions such that the first point is positive
    F = X.T[:g_naxes]
    for i in range(g_naxes):
        if F[i][0] < 0:
            F[i] *= -1
    # create the image string
    image_string = create_image_string(fs.imageformat, physical_size, F, fs.xaxis_length)
    # create the response headers
    image_basename = 'mds'
    response_headers = []
    format_to_content_type = {'svg':'image/svg+xml', 'png':'image/png', 'pdf':'application/pdf', 'ps':'application/postscript'}
    response_headers.append(('Content-Type', format_to_content_type[fs.imageformat]))
    response_headers.append(('Content-Disposition', "%s; filename=%s.%s" % (fs.contentdisposition, image_basename, fs.imageformat)))
    # return the response
    return response_headers, image_string


def create_image_string(image_format, physical_size, F, xaxis_length):
    """
    This function is about drawing the tree.
    param image_format: the image extension
    @param physical_size: the width and height of the image in pixels
    @param F: eigenfunction samples
    @param xaxis_length: True if plotting vs length; False if plotting vs PC1
    @return: the image as a string
    """
    # before we begin drawing we need to create the cairo surface and context
    cairo_helper = CairoUtil.CairoHelper(image_format)
    surface = cairo_helper.create_surface(physical_size[0], physical_size[1])
    context = cairo.Context(surface)
    # define some helper variables
    x0 = physical_size[0] / 2.0
    y0 = physical_size[1] / 2.0
    # draw an off-white background
    context.save()
    context.set_source_rgb(.9, .9, .9)
    context.paint()
    context.restore()
    # draw the axes which are always in the center of the image
    context.save()
    context.set_source_rgb(0, 0, 0)
    context.move_to(x0, 0)
    context.line_to(x0, physical_size[1])
    context.stroke()
    context.move_to(0, y0)
    context.line_to(physical_size[0], y0)
    context.stroke()
    context.restore()
    # define red, blue, green, orange
    colors = [
            (1.0, 0.0, 0.0),
            (0.0, 0.0, 1.0),
            (0.0, 1.0, 0.0),
            (1.0, 0.5, 0.0)]
    # add the sum of the eigenfunctions with a gray color
    #F = np.vstack([F, np.sum(F, 0)])
    #colors.append((0.5, 0.5, 0.5))
    # draw the eigenfunctions
    for color, v in zip(colors, F):
        # define the sequence of physical points
        seq = []
        for i, y in enumerate(v):
            xprogress = i / (len(v) - 1.0)
            if xaxis_length:
                x = x0 + (xprogress - 0.5) * 0.75 * physical_size[0]
            else:
                x = x0 + F[0][i] * physical_size[0]
            y = y0 + y * physical_size[1]
            seq.append((x, y))
        # draw the eigenfunction
        context.save()
        context.set_source_rgb(*color)
        for (xa, ya), (xb, yb) in iterutils.pairwise(seq):
            context.move_to(xa, ya)
            context.line_to(xb, yb)
            context.stroke()
        context.restore()
    # create the image
    return cairo_helper.get_image_string()

