"""Plot some connected points.
"""

from StringIO import StringIO
import math

import cairo

from SnippetUtil import HandlingError
import Codon
import Form
import FormOut
import CairoUtil
import Util
import const

g_default_data = const.read('20100730f')

def get_form():
    """
    @return: the body of a form
    """
    # define the default labeled points
    lines = Util.get_stripped_lines(g_default_data.splitlines())
    # define the form objects
    form_objects = [
            Form.MultiLine('graph_data', 'connected points',
                '\n'.join(lines)),
            Form.Integer('total_width', 'total image width',
                640, low=3, high=2000),
            Form.Integer('total_height', 'total image height',
                480, low=3, high=2000),
            Form.Integer('border', 'image border size',
                10, low=0, high=2000),
            Form.ImageFormat()]
    return form_objects

def get_form_out():
    return FormOut.Image('plot')

def get_image_string(points, edges, t_width, t_height, border, image_format):
    """
    @param points: an ordered list of (x, y) pairs
    @param edges: a set of point index pairs
    @param t_width: the image width in pixels
    @param t_height: the image height in pixels
    @param border: the width and height of the image border in pixels
    @param image_format: the requested format of the image
    """
    width = t_width - 2*border
    height = t_height - 2*border
    assert width >= 1
    assert height >= 1
    # get the x and y coordinates of the points
    x_coords, y_coords_raw = zip(*points)
    # Flip the y coordinates so that greater values of y
    # are shown near the top of the image.
    y_coords = [-y for y in y_coords_raw]
    unscaled_width = max(x_coords) - min(x_coords)
    unscaled_height = max(y_coords) - min(y_coords)
    # rescale the points to fit in the drawable part of the image
    c_width = width / unscaled_width
    c_height = height / unscaled_height
    c = min(c_width, c_height)
    x_rescaled = [x*c for x in x_coords]
    y_rescaled = [y*c for y in y_coords]
    # translate the points so that their minimum coordinate is zero
    x_min = min(x_rescaled)
    y_min = min(y_rescaled)
    x_trans = [x-x_min for x in x_rescaled]
    y_trans = [y-y_min for y in y_rescaled]
    # translate the points to account for the border
    x_final = [x+border for x in x_trans]
    y_final = [y+border for y in y_trans]
    # create the surface
    cairo_helper = CairoUtil.CairoHelper(image_format)
    surface = cairo_helper.create_surface(t_width, t_height)
    context = cairo.Context(surface)
    # draw the background
    context.save()
    context.set_source_rgb(.9, .9, .9)
    context.paint()
    context.restore()
    # draw the points
    radius = 2.0
    for x, y in zip(x_final, y_final):
        # draw a filled circle
        context.save()
        context.arc(x, y, radius, 0, 2 * math.pi)
        context.close_path()
        context.fill()
        context.restore()
    # draw the edges
    for i, j in edges:
        context.save()
        context.move_to(x_final[i], y_final[i])
        context.line_to(x_final[j], y_final[j])
        context.close_path()
        context.stroke()
        context.restore()
    # get the image string
    return cairo_helper.get_image_string()

def read_points_and_edges(multiline):
    """
    @param multiline: input like the default data
    @return: a list of (x, y) points and a set of point-index-pair edges
    """
    lines = Util.get_stripped_lines(multiline.splitlines())
    try:
        POINTS_index = lines.index('POINTS')
    except ValueError as e:
        raise HandlingError('expected a line that says POINTS')
    try:
        EDGES_index = lines.index('EDGES')
    except ValueError as e:
        raise HandlingError('expected a line that says EDGES')
    if POINTS_index > EDGES_index:
        raise HandlingError('expected points before edges')
    # read the points
    points = []
    points_lines = lines[POINTS_index + 1: EDGES_index]
    for i, point_line in enumerate(points_lines):
        values = point_line.split()
        if len(values) != 3:
            raise ValueError('each POINTS row should have three values')
        s_index, s_x, s_y = values
        try:
            index = int(s_index)
            x = float(s_x)
            y = float(s_y)
        except ValueError as e:
            raise HandlingError('a value in a POINTS row has the wrong type')
        if index != i:
            raise HandlingError('the POINTS indices should match their order')
        point = (x, y)
        points.append(point)
    # read the edges
    edges = set()
    edges_lines = lines[EDGES_index + 1:]
    for edge_line in edges_lines:
        values = edge_line.split()
        if len(values) != 2:
            raise ValueError('each EDGES row should have three values')
        s_i, s_j = values
        try:
            i = int(s_i)
            j = int(s_j)
        except ValueError as e:
            raise HandlingError('a value in a EDGES row has the wrong type')
        edge = (i, j)
        edges.add(edge)
    # return the points and edges
    return points, edges


def get_response_content(fs):
    # read the points and edges
    points, edges = read_points_and_edges(fs.graph_data)
    # get the width and height of the drawable area of the image
    width = fs.total_width - 2*fs.border
    height = fs.total_height - 2*fs.border
    if width < 1 or height < 1:
        msg = 'the image dimensions do not allow for enough drawable area'
        raise HandlingError(msg)
    # draw the image
    try:
        ext = Form.g_imageformat_to_ext[fs.imageformat]
        return get_image_string(points, edges,
                fs.total_width, fs.total_height, fs.border, ext)
    except CairoUtil.CairoUtilError as e:
        raise HandlingError(e)
