"""Plot points using resistance distance MDS of degree 1 nodes.

The image is constrained to have the center at the origin.
Axes may be added.
Points may be labeled according by index.
"""

from StringIO import StringIO
import math
import os

import cairo
import numpy as np

from SnippetUtil import HandlingError
import Form
import FormOut
import CairoUtil
import Euclid
import BuildTreeTopology
import SchurAlgebra
import Util
import const

g_data = const.read('20100125a')

def get_default_data():
    lines = Util.get_stripped_lines(g_data.splitlines())
    return '\n'.join(lines)

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
    def __init__(self, width, height,
            all_black, show_labels,
            axis_info, border_info, image_format):
        self.width = width
        self.height = height
        self.all_black = all_black
        self.show_labels = show_labels
        self.axis_info = axis_info
        self.border_info = border_info
        self.image_format = image_format

def get_image_string(x_coords, y_coords,
        point_colors, edges, image_info, scaling_factor):
    """
    @param edges: ignore this for now
    """
    # unpack the total width and height of the image
    t_width = image_info.width
    t_height = image_info.height
    # unpack axis options
    show_x = image_info.axis_info.show_x
    show_y = image_info.axis_info.show_y
    show_labels = image_info.show_labels
    all_black = image_info.all_black
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
    dot_radius = 3.0
    for x, y, point_color in zip(x_coords, y_coords, point_colors):
        context.save()
        if not all_black:
            context.set_source_rgb(*point_color)
        context.arc(x, y, dot_radius, 0, 2 * math.pi)
        context.close_path()
        context.fill()
        context.restore()
    # Draw the labels.
    if show_labels:
        labels = [str(i) for i, x in enumerate(x_coords)]
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
            Form.MultiLine('graph_data', 'point data', get_default_data()),
            Form.CheckGroup('axis_options', 'axis options', [
                Form.CheckItem('flip_x', 'flip x axis'),
                Form.CheckItem('flip_y', 'flip y axis'),
                Form.CheckItem('show_x', 'show x axis', True),
                Form.CheckItem('show_y', 'show y axis', True)]),
            Form.CheckGroup('vis_options', 'visualization options', [
                Form.CheckItem('show_labels', 'show labels')]),
            Form.RadioGroup('edge_weight_options', 'edge weights', [
                Form.RadioItem('unweighted', 'all weights are 1.0', True),
                Form.RadioItem('weighted', 'weights are inverse distances')]),
            Form.RadioGroup('color_options', 'color options', [
                Form.RadioItem('black', 'all black'),
                Form.RadioItem('color', 'first axis coloration', True)]),
            Form.CheckGroup('schur_options', 'Schur options', [
                Form.CheckItem('schur', 'use schur complement', True)]),
            Form.Float('scale', 'scaling factor',
                '150.0', low_exclusive=0),
            Form.Integer('width', 'image width in pixels',
                640, low=1, high=2000),
            Form.Integer('height', 'image height in pixels',
                480, low=1, high=2000),
            Form.ImageFormat(),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.Image('plot')

def edges_to_node_degrees(edges):
    """
    @param edges: a set of point-index-pair edges
    @return: a list of node degrees
    """
    edge_list = list(edges)
    first_indices, second_indices = zip(*edge_list)
    index_set = set(first_indices) | set(second_indices)
    npoints = max(index_set) + 1
    index_list = list(sorted(index_set))
    if index_list != range(npoints):
        raise ValueError('degree zero nodes are not allowed')
    index_to_count = [0]*npoints
    for i, j in edges:
        index_to_count[i] += 1
        index_to_count[j] += 1
    return index_to_count

def read_points_and_edges(multiline):
    """
    @param multiline: input like the default data
    @return: a list of (x, y) points and a set of point-index-pair edges
    """
    lines = [x.strip() for x in StringIO(multiline).readlines()]
    lines = [x for x in lines if x]
    try:
        POINTS_index = lines.index('POINTS')
    except ValueError as e
        raise HandlingError('expected a line that says POINTS')
    try:
        EDGES_index = lines.index('EDGES')
    except ValueError as e
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
        except ValueError as e
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
        except ValueError as e
            raise HandlingError('a value in a EDGES row has the wrong type')
        edge = (i, j)
        edges.add(edge)
    # return the points and edges
    return points, edges

def valuations_to_colors(valuations):
    """
    @param valuations: Principal axis valuations of the points
    @return: a conformant list of cairo formatted rgb triples
    """
    dark_blue = (0.0, 0.0, 1.0)
    light_blue = (0.4, 0.4, 1.0)
    dark_red = (1.0, 0.0, 0.0)
    light_red = (1.0, 0.4, 0.4)
    negative_valuations = [v for v in valuations if v <= 0]
    positive_valuations = [v for v in valuations if v > 0]
    neg_median = np.median(negative_valuations)
    pos_median = np.median(positive_valuations)
    colors = []
    for v in valuations:
        if v <= 0:
            if v <= neg_median:
                c = dark_blue
            else:
                c = light_blue
        else:
            if v >= pos_median:
                c = dark_red
            else:
                c = light_red
        colors.append(c)
    return colors

def edges_to_laplacian(edges, edge_weights):
    """
    Return a Laplacian matrix given a list of unweighted edges.
    Vertices of the graph are expected to be indexed
    sequentially starting at zero.
    @param edges: a list of connected point index pairs
    @param edge_weights: a list of edge weights
    @return: a numpy array laplacian matrix
    """
    source_indices, sink_indices = zip(*edges)
    all_indices = set(source_indices) | set(sink_indices)
    npoints = max(all_indices) + 1
    if all_indices != set(range(npoints)):
        raise ValueError('the graph is not connected')
    # create the adjacency matrix
    A = np.zeros((npoints, npoints))
    for weight, (i, j) in zip(edge_weights, edges):
        A[i, j] = weight
        A[j, i] = weight
    # create the laplacian matrix
    L = Euclid.adjacency_to_laplacian(A)
    return L

def get_response_content(fs):
    # use a default border; actually this is ignored
    border_info = BorderInfo(10, 10)
    # Collect the image format information.
    axis_info = AxisInfo(fs.flip_x, fs.flip_y, fs.show_x, fs.show_y)
    # read the points and edges
    points, edges = read_points_and_edges(fs.graph_data)
    # define edge weights
    if fs.weighted:
        np_points = [np.array(p) for p in points]
        dists = [np.linalg.norm(np_points[j] - np_points[i]) for i, j in edges]
        weights = [1.0 / d for d in dists]
    else:
        weights = [1.0 for e in edges]
    # create the full laplacian
    L = edges_to_laplacian(edges, weights)
    if fs.schur:
        # remove internal nodes by schur complementation
        index_to_degree = edges_to_node_degrees(edges)
        internal_indices = set(i
                for i, d in enumerate(index_to_degree) if d > 1)
        L_schur = SchurAlgebra.mschur(L, internal_indices)
        L_final = L_schur
    else:
        L_final = L
    # define the point colors using the graph Fiedler loadings
    G = np.linalg.pinv(L_final)
    X = Euclid.dccov_to_points(G)
    points = [(p[0], p[1]) for p in X]
    xs, ys = zip(*points)
    colors = valuations_to_colors(xs)
    # draw the image
    ext = Form.g_imageformat_to_ext[fs.imageformat]
    image_info = ImageInfo(fs.width, fs.height,
            fs.black, fs.show_labels,
            axis_info, border_info, ext)
    return get_image_string(xs, ys, colors, edges, image_info, fs.scale)
