"""Plot labeled points with origin at center and preserving aspect ratio.

For example, this can be used for a city MDS plot.
"""

import StringIO
import math

import cairo

from SnippetUtil import HandlingError
import Form
import CairoUtil

g_default_data = """Boston  -905.527460092  666.956097193
Seattle 1573.80442188   425.46464836
Phoenix 1030.45352721   -564.133058397
Miami   -948.765620649  -600.926346854
Raleigh -749.964868349  72.6386596982"""


class BorderInfo:
    def __init__(self, border_x, border_y):
        self.border_x = border_x
        self.border_y = border_y

class AxisInfo:
    def __init__(self, flip_x, flip_y, show_x, show_y):
        self.flip_x = flip_x
        self.flip_y = flip_y
        self.show_x = show_x
        self.show_y = show_y

class ImageInfo:
    def __init__(self, width, height, axis_info, border_info, image_format):
        self.width = width
        self.height = height
        self.axis_info = axis_info
        self.border_info = border_info
        self.image_format = image_format

def get_image_string(x_coords, y_coords, labels, image_info):
    # unpack the total width and height of the image
    t_width = image_info.width
    t_height = image_info.height
    # unpack axis visualization options
    show_x = image_info.axis_info.show_x
    show_y = image_info.axis_info.show_y
    # flip axes if necessary
    if image_info.axis_info.flip_x:
        x_coords = [-x for x in x_coords]
    if image_info.axis_info.flip_y:
        y_coords = [-y for y in y_coords]
    # unpack border info
    border_x = image_info.border_info.border_x
    border_y = image_info.border_info.border_y
    # unpack the image format
    image_format = image_info.image_format
    # Define the scaling factors in the x+, x-, y+, and y- directions
    # which would fill the entire drawable (non-border) region.
    # The smallest of these scaling factors will be used for all directions.
    xpos_coords = [x for x in x_coords if x > 0]
    xneg_coords = [x for x in x_coords if x < 0]
    ypos_coords = [y for y in y_coords if y > 0]
    yneg_coords = [y for y in y_coords if y < 0]
    sf_list = []
    if xpos_coords:
        available = (t_width / 2.0) - border_x
        used = max(xpos_coords)
        sf_list.append(available / used)
    if xneg_coords:
        available = (t_width / 2.0) - border_x
        used = -min(xneg_coords)
        sf_list.append(available / used)
    if ypos_coords:
        available = (t_height / 2.0) - border_y
        used = max(ypos_coords)
        sf_list.append(available / used)
    if yneg_coords:
        available = (t_height / 2.0) - border_y
        used = -min(yneg_coords)
        sf_list.append(available / used)
    scaling_factor = min(sf_list)
    # Update the coordinates using the new scaling factor.
    x_coords = [x*scaling_factor for x in x_coords]
    y_coords = [y*scaling_factor for y in y_coords]
    # Translate the coordinates
    # so that the origin is at the center of the image.
    x_coords = [x + t_width / 2.0 for x in x_coords]
    y_coords = [y + t_height / 2.0 for y in y_coords]
    # create the surface
    cairo_helper = CairoUtil.CairoHelper(image_format)
    surface = cairo_helper.create_surface(t_width, t_height)
    context = cairo.Context(surface)
    # draw the background
    context.save()
    context.set_source_rgb(.9, .9, .9)
    context.paint()
    context.restore()
    # Draw the x axis if requested.
    if show_x:
        context.save()
        context.set_source_rgb(.9, .7, .7)
        context.move_to(0, t_height / 2.0)
        context.line_to(t_width, t_height / 2.0)
        context.stroke()
        context.restore()
    # Draw the y axis if requested.
    if show_y:
        context.save()
        context.set_source_rgb(.9, .7, .7)
        context.move_to(t_width / 2.0, 0)
        context.line_to(t_width / 2.0, t_height)
        context.stroke()
        context.restore()
    # Draw the points.
    dot_radius = 2.0
    for x, y in zip(x_coords, y_coords):
        context.save()
        context.arc(x, y, dot_radius, 0, 2 * math.pi)
        context.close_path()
        context.fill()
        context.restore()
    # Draw the labels.
    for label, x, y in zip(labels, x_coords, y_coords):
        context.save()
        context.move_to(x, y)
        context.show_text(label)
        context.restore()
    # get the image string
    return cairo_helper.get_image_string()

def get_form():
    """
    @return: a list of form objects
    """
    # define the form objects
    form_objects = [
            Form.MultiLine('points', 'labeled points', g_default_data),
            Form.CheckGroup('axis_options', 'axis options', [
                Form.CheckItem('flip_x', 'flip x axis'),
                Form.CheckItem('flip_y', 'flip y axis'),
                Form.CheckItem('show_x', 'show x axis'),
                Form.CheckItem('show_y', 'show y axis')]),
            Form.Integer('border_x', 'padding on the left and right',
                10, low=0),
            Form.Integer('border_y', 'padding on the top and bottom',
                10, low=0),
            Form.Integer('width', 'image width in pixels',
                640, low=1, high=2000),
            Form.Integer('height', 'image height in pixels',
                480, low=1, high=2000),
            Form.RadioGroup('imageformat', 'image format', [
                Form.RadioItem('png', 'png', True),
                Form.RadioItem('svg', 'svg'),
                Form.RadioItem('pdf', 'pdf'),
                Form.RadioItem('ps', 'ps')]),
            Form.RadioGroup('contentdisposition', 'image delivery', [
                Form.RadioItem('inline', 'view the image', True),
                Form.RadioItem('attachment', 'download the image')])]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # Collect the image format information.
    border_info = BorderInfo(fs.border_x, fs.border_y)
    axis_info = AxisInfo(fs.flip_x, fs.flip_y, fs.show_x, fs.show_y)
    w, h = fs.width, fs.height
    image_info = ImageInfo(w, h, axis_info, border_info, fs.imageformat)
    # Parse the label and point input.
    lines = [x.strip() for x in StringIO.StringIO(fs.points).readlines()]
    lines = [x for x in lines if x]
    triples = [x.split() for x in lines]
    parsed_triples = []
    for triple in triples:
        if len(triple) != 3:
            raise HandlingError('expected three entries per row')
        label, sx, sy = triple
        try:
            x = float(sx)
            y = float(sy)
        except ValueError, e:
            raise HandlingError('expected x and y coordinates to be numbers')
        parsed_triples.append((label, x, y))
    labels, X, Y = zip(*parsed_triples)
    # Get the image.
    image_string = get_image_string(X, Y, labels, image_info)
    # start writing the response type
    response_headers = []
    # specify the content type
    format_to_content_type = {
            'svg':'image/svg+xml', 'png':'image/png',
            'pdf':'application/pdf', 'ps':'application/postscript'}
    content_type = format_to_content_type[fs.imageformat]
    response_headers.append(('Content-Type', content_type))
    # specify the content disposition
    image_filename = 'plot.' + fs.imageformat
    content_disp = '%s; filename=%s' % (fs.contentdisposition, image_filename)
    header = ('Content-Disposition', content_disp)
    response_headers.append(header)
    # return the response
    return response_headers, image_string
