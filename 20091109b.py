"""Plot the first few eigenvectors of the Laplacian of a path.

The eigenvectors are plotted in a way
that makes them look like eigenfunction approximations.
"""

import numpy as np
import cairo

from SnippetUtil import HandlingError
import Form
import FormOut
import Euclid
import CairoUtil
import iterutils

def get_form():
    """
    @return: a list of form objects
    """
    # define the form objects
    form_objects = [
            Form.Integer('nvertices', 'use this many vertices',
                20, low=2, high=100),
            Form.Integer('naxes', 'plot this many eigenvectors',
                4, low=1),
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
            Form.ImageFormat(),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.Image('eigenfunctions')

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

def get_response_content(fs):
    # check input compatibility
    if fs.nvertices < fs.naxes+1:
        msg_a = 'attempting to plot too many eigenvectors '
        msg_b = 'for the given number of vertices'
        raise ValueError(msg_a + msg_b)
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
    F = X.T[:fs.naxes]
    for i in range(fs.naxes):
        if F[i][0] < 0:
            F[i] *= -1
    # draw the image
    try:
        ext = Form.g_imageformat_to_ext[fs.imageformat]
        return create_image_string(ext, physical_size, F, fs.xaxis_length)
    except CairoUtil.CairoUtilError as e
        raise HandlingError(e)

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
    for i, v in enumerate(F):
        # define the color
        if i < len(colors):
            color = colors[i]
        else:
            color = (0,0,0)
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

