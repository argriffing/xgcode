"""Plot some labeled points.
"""

from StringIO import StringIO
import math

import cairo

from SnippetUtil import HandlingError
import Codon
import Form
import Util
import CairoUtil

def todec(degrees, minutes):
    """
    @param degrees: latitude or longitude degrees
    @param minutes: latitude or longitude minutes
    @return: a floating point number
    """
    return degrees + minutes / 60.0

def get_form():
    """
    @return: the body of a form
    """
    # define the default labeled points
    default_labeled_points = [
            ('raleigh', -todec(78, 39), todec(35, 46)),
            ('atl', -todec(84, 23), todec(33, 45)),
            ('dc', -todec(77, 02), todec(38, 53)),
            ('pitt', -todec(79, 57), todec(40, 27)),
            ('philly', -todec(75, 10), todec(39, 57))]
    tsv = ['\t'.join([(label, str(x), str(y)])
        for label, x, y in default_labeled_points]
    default_labeled_points_string = '\n'.join(tsv)
    # define the form objects
    form_objects = [
            Form.MultiLine('labeled_points', 'labeled points',
                default_labeled_points_string),
            Form.Integer('total_width', 'total image width',
                640, low=3, high=2000),
            Form.Integer('total_height', 'total image height',
                480, low=3, high=2000),
            Form.Integer('border', 'image border size',
                10, low=0, high=2000),
            Form.RadioGroup('imageformat', 'image format options', [
                Form.RadioItem('png', 'png', True),
                Form.RadioItem('svg', 'svg'),
                Form.RadioItem('pdf', 'pdf'),
                Form.RadioItem('ps', 'ps')]),
            Form.RadioGroup('contentdisposition', 'image delivery options', [
                Form.RadioItem('inline', 'view the image', True),
                Form.RadioItem('attachment', 'download the image')])]
    return form_objects

def get_image_string(labels, points, total_width, total_height, border, image_format):
    """
    @param labels: an ordered list of point labels
    @param points: an ordered list of (x, y) points
    @param total_width: the image width in pixels
    @param total_height: the image height in pixels
    @param border: the width and height of the image border in pixels
    @param image_format: the requested format of the image
    """
    assert len(labels) == len(points)
    width = total_width - 2*border
    height = total_height - 2*border
    assert width >= 1
    assert height >= 1
    # get the x and y coordinates of the points
    x_coords, y_coords_raw = zip(*points)
    # flip the y coordinates so that greater values of y are shown near the top of the image
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
    surface = cairo_helper.create_surface(total_width, total_height)
    context = cairo.Context(surface)
    # draw the background
    context.save()
    context.set_source_rgb(.9, .9, .9)
    context.paint()
    context.restore()
    # draw the labeled points
    radius = 2.0
    for label, x, y in zip(labels, x_final, y_final):
        # draw a filled circle
        context.save()
        context.arc(x, y, radius, 0, 2 * math.pi)
        context.close_path()
        context.fill()
        context.restore()
        # draw the label
        context.move_to(x, y)
        context.show_text(label)
        context.close_path()
        context.stroke()
    # get the image string
    return cairo_helper.get_image_string()

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # read the labeled points
    labels = []
    points = []
    labeled_point_lines = Util.get_stripped_lines(StringIO(fs.labeled_points))
    for line in labeled_point_lines:
        labeled_point = line.split()
        if len(labeled_point) != 3:
            raise HandlingError('each line should have three whitespace separated elements')
        label, x_string, y_string = labeled_point
        try:
            x = float(x_string)
            y = float(y_string)
        except ValueError, e:
            raise HandlingError('expected the coordinates to be floating point numbers')
        labels.append(label)
        points.append((x, y))
    # get the width and height of the drawable area of the image
    width = fs.total_width - 2*fs.border
    height = fs.total_height - 2*fs.border
    if width < 1 or height < 1:
        raise HandlingError('the image dimensions do not allow for enough drawable area')
    # draw the image
    try:
        image_string = get_image_string(labels, points, fs.total_width, fs.total_height, fs.border, fs.imageformat)
    except CairoUtil.CairoUtilError, e:
        raise HandlingError(e)
    # begin the response
    response_headers = []
    # specify the content type
    format_to_content_type = {'svg':'image/svg+xml', 'png':'image/png', 'pdf':'application/pdf', 'ps':'application/postscript'}
    response_headers.append(('Content-Type', format_to_content_type[fs.imageformat]))
    # specify the content disposition
    image_filename = 'plot.' + fs.imageformat
    response_headers.append(('Content-Disposition', "%s; filename=%s" % (fs.contentdisposition, image_filename)))
    # return the response
    return response_headers, image_string
