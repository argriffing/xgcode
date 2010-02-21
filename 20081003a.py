"""Demonstrate kernel PCA.

This is not so useful.
"""

from StringIO import StringIO
import math

import numpy as np
import cairo

from SnippetUtil import HandlingError
import Form
import CairoUtil
import MatrixUtil
import Clustering
import SpiralSampler

def get_form():
    """
    @return: a list of form objects
    """
    # sample some points to use as the defaults
    M = np.array(list(SpiralSampler.gen_points(200, .01)))
    # define the form objects
    form_objects = [
            Form.Matrix('points', 'points', M),
            Form.RadioGroup('imageformat', 'image format', [
                Form.RadioItem('png', 'png', True),
                Form.RadioItem('svg', 'svg'),
                Form.RadioItem('pdf', 'pdf'),
                Form.RadioItem('ps', 'ps')]),
            Form.RadioGroup('contentdisposition', 'image delivery', [
                Form.RadioItem('inline', 'view the image', True),
                Form.RadioItem('attachment', 'download the image')])]
    return form_objects

def get_eigenvectors(row_major_matrix):
    """
    This gets a couple of left eigenvectors
    because of the standard format of rate matrices.
    @param row_major_matrix: this is supposed to be a rate matrix
    @return: a pair of eigenvectors
    """
    R = np.array(row_major_matrix)
    w, vl, vr = np.linalg.eig(R, left=True, right=True)
    eigenvalue_info = list(sorted((abs(x), i) for i, x in enumerate(w)))
    stationary_eigenvector_index = eigenvalue_info[0][1]
    first_axis_eigenvector_index = eigenvalue_info[1][1]
    second_axis_eigenvector_index = eigenvalue_info[2][1]
    va = vl.T[first_axis_eigenvector_index]
    vb = vl.T[second_axis_eigenvector_index]
    return va, vb

def get_rescaled_vector(v):
    """
    @param v: an array or list of floating point values
    @return: a list of floating point values rescaled to be in the range (0, 1)
    """
    width = max(v) - min(v)
    if not width:
        return [.5 for x in v]
    return [(x-min(v)) / width for x in v]

def get_image_helper(xmin, ymin, xmax, ymax, vx, vy, image_size, image_format):
    """
    Get an image string.
    Input comprises excruciating details about the image
    and its scaling properties.
    @param xmin: the smallest x value allowed in the viewport
    @param ymin: the smallest y value allowed in the viewport
    @param xmax: the greatest x value allowed in the viewport
    @param ymax: the greatest y value allowed in the viewport
    @param vx: the list of x coordinates of the points
    @param vy: the list of y coordinates of the points
    @param image_size: the width and height of the image in pixels
    @param image_format: like 'svg', 'png', 'ps', 'pdf', et cetera
    @return: an image string
    """
    width, height = image_size
    n = len(vx)
    # rescale the x and y vectors to be between zero and one
    xextent = xmax - xmin
    yextent = ymax - ymin
    rescaled_vx = [(x - xmin) / xextent for x in vx]
    rescaled_vy = [(y - ymin) / yextent for y in vy]
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
    # begin drawing the axes
    context.save()
    context.set_source_rgb(.9, .7, .7)
    # draw the y axis
    tx = -xmin/xextent
    xzero = (tx * (1 - 2*border_fraction) + border_fraction) * width
    context.move_to(xzero, 0)
    context.line_to(xzero, height)
    context.stroke()
    # draw the x axis
    ty = -ymin/yextent
    yzero = (ty * (1 - 2*border_fraction) + border_fraction) * height
    context.move_to(0, yzero)
    context.line_to(width, yzero)
    context.stroke()
    # stop drawing the axes
    context.restore()
    # draw a scatter plot of the states using the eigenvectors as axes
    for i, (x, y) in enumerate(zip(rescaled_vx, rescaled_vy)):
        if i < n/2:
            state_string = 'x'
        else:
            state_string = 'o'
        nx = (x * (1 - 2*border_fraction) + border_fraction) * width
        ny = (y * (1 - 2*border_fraction) + border_fraction) * height
        context.move_to(nx, ny)
        context.show_text(state_string)
    # get the image string
    return cairo_helper.get_image_string()

def get_image(row_major_matrix, width_and_height, image_format, draw_axes):
    """
    @param row_major_matrix: this is supposed to be a rate matrix
    @param width_and_height: the dimensions of the output image
    @param image_format: like 'svg', 'png', 'ps', 'pdf', et cetera
    @param draw_axes: True if axes should be drawn
    @return: a string containing the image data
    """
    # FIXME axes are always drawn
    va, vb = get_eigenvectors(row_major_matrix)
    xmin, xmax = min(va), max(va)
    ymin, ymax = min(vb), max(vb)
    return get_image_helper(xmin, ymin, xmax, ymax,
            va, vb, width_and_height, image_format)

def points_to_distance_matrix(points):
    """
    @param points: an numpy array of points in euclidean space
    @return: a distance matrix
    """
    n = len(points)
    D = np.zeros((n,n))
    for i, point_a in enumerate(points):
        for j, point_b in enumerate(points):
            v = point_b - point_a
            D[i][j] = np.dot(v, v)
    return D

def make_laplacian(D, stddev):
    """
    @param D: a distance matrix
    @param stddev: the standard deviation of the gaussian kernel
    """
    n = len(D)
    # get the similarity matrix
    similarity_matrix = []
    for i in range(n):
        row = []
        for j in range(n):
            if i == j:
                similarity = 0
            else:
                distance = D[i][j]
                similarity = (1/stddev)*math.exp(-(distance / stddev)**2)
            row.append(similarity)
        similarity_matrix.append(row)
    # get the laplacian matrix
    laplacian_matrix = [row[:] for row in similarity_matrix]
    for i in range(n):
        laplacian_matrix[i][i] -= sum(laplacian_matrix[i])
    return laplacian_matrix

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # get the points
    points = fs.points
    # get the distance matrix from the euclidean points
    n = len(points)
    D = points_to_distance_matrix(points)
    # get the similarity matrix
    #stddev = 0.2 # this gives a good separation
    #stddev = 0.8 # this gives a less good separation
    #stddev = 1.0 # this gives something that looks like the original matrix
    stddev = 0.3
    laplacian_matrix = make_laplacian(D, stddev)
    # create the file name
    image_format = fs.imageformat
    # draw the image
    try:
        image_size = (640, 480)
        draw_axes = True
        image_string = get_image(laplacian_matrix, image_size,
                image_format, draw_axes)
    except CairoUtil.CairoUtilError, e:
        raise HandlingError(e)
    # start writing the response type
    response_headers = []
    # specify the content type
    format_to_content_type = {'svg':'image/svg+xml', 'png':'image/png', 'pdf':'application/pdf', 'ps':'application/postscript'}
    response_headers.append(('Content-Type', format_to_content_type[image_format]))
    # specify the content disposition
    image_extension = image_format
    image_filename = 'scatterplot' + '.' + image_extension
    response_headers.append(('Content-Disposition', "%s; filename=%s" % (fs.contentdisposition, image_filename)))
    # return the response
    return response_headers, image_string

def main():
    """
    Make images from the command line.
    """
    # set the sampling options
    npoints_per_group = 200
    sample_stddev = .05
    # video options
    nframes = 51
    # clustering options
    kernel_stddev_initial = 1.0
    kernel_stddev_final = 0.1
    kernel_stddev_gap = kernel_stddev_final - kernel_stddev_initial
    kernel_stddev_step = kernel_stddev_gap / (nframes - 1)
    # sample some points
    points = np.array(list(SpiralSampler.gen_points(
        npoints_per_group, sample_stddev)))
    # make the distance matrix
    D = points_to_distance_matrix(points)
    # get two eigenvectors for each frame
    frame_eigenvectors = []
    for i in range(nframes):
        print 'calculating eigendecomposition for frame', i
        # set the kernel standard deviation for this frame
        kernel_stddev = kernel_stddev_initial + i * kernel_stddev_step
        # make the laplacian matrix
        laplacian_matrix = make_laplacian(D, kernel_stddev)
        # get the eigenvectors
        eigenvectors = get_eigenvectors(laplacian_matrix)
        # possibly flip the signs of the eigenvectors
        if i:
            vx, vy = eigenvectors
            last_vx, last_vy = frame_eigenvectors[i-1]
            # possibly flip the x vector
            match_count = sum(1 for a, b in zip(vx, last_vx) if a*b >= 0)
            if match_count < len(points) / 2:
                print 'flipping x'
                vx = [-value for value in vx]
            # possibly flip the y vector
            match_count = sum(1 for a, b in zip(vy, last_vy) if a*b >= 0)
            if match_count < len(points) / 2:
                print 'flipping y'
                vy = [-value for value in vy]
            eigenvectors = (vx, vy)
        frame_eigenvectors.append(eigenvectors)
    # get all the eigenvector elements in each direction
    x_eigenvector_elements = []
    y_eigenvector_elements = []
    for x_eigenvector, y_eigenvector in frame_eigenvectors:
        x_eigenvector_elements.extend(x_eigenvector)
        y_eigenvector_elements.extend(y_eigenvector)
    # get the max and min of the eigenvectors in each direction
    xmin, xmax = min(x_eigenvector_elements), max(x_eigenvector_elements)
    ymin, ymax = min(y_eigenvector_elements), max(y_eigenvector_elements)
    # write the files
    for i, (vx, vy) in enumerate(frame_eigenvectors):
        print 'writing frame', i
        # get the image
        image_format = 'png'
        image_size = (800, 600)
        draw_axes = True
        image_string = get_image_helper(xmin, ymin, xmax, ymax,
                vx, vy, image_size, image_format)
        # get the filename
        filename = 'frame-%04d.png' % i
        # write the image
        fout = open(filename, 'wb')
        fout.write(image_string)
        fout.close()


if __name__ == '__main__':
    main()
