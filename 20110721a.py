"""
Draw a tikz illustration of harmonic interpolation.
"""

import numpy as np

import math

import Form
import FormOut
import tikz

def get_form():
    """
    @return: a list of form objects
    """
    # define the form objects
    form_objects = [
            Form.Float('d1', 'd1', 1.0),
            Form.Float('d2', 'd2', 2.0),
            Form.Float('d3', 'd3', 3.0),
            Form.Float('y1', 'y1', 1.0),
            Form.Float('y2', 'y2', 2.5),
            Form.Float('y3', 'y3', 3.5),
            Form.TikzFormat()]
    return form_objects

def get_form_out():
    return FormOut.Tikz()

def rotate_to_view(p):
    """
    Rotate a few degrees around the z axis and then around the new y axis.
    @param p: a 3d point
    @return: a 3d point rotated around the origin
    """
    # use pi/4 for a more standard rotation
    theta = -math.pi / 12
    c = math.cos(theta)
    s = math.sin(theta)
    x0, y0, z0 = p
    # first rotation
    x1 = x0 * c - y0 * s
    y1 = y0 * c + x0 * s
    z1 = z0
    # second rotation
    x2 = x1 * c - z1 * s
    y2 = y1
    z2 = z1 * c + x1 * s
    return np.array([x2, y2, z2])

def get_domain_vertex_tikz_lines(points):
    """
    Define the domain vertices.
    @param points: a sequence of four 2D numpy points representing the domain
    """
    arr = []
    # define the dot style
    dotstyle = 'fill=black,circle,inner sep=0pt,minimum size=3pt'
    arr.append('\\tikzstyle{dot}=[%s]' % dotstyle)
    # draw the dots
    for i, (x, y) in enumerate(points):
        style = ''
        name = 'n%s' % i
        line = '\\node[dot] (%s) at (%04f, %04f) {};' % (name, x, y)
        arr.append(line)
    return arr

def get_domain_edge_tikz_lines():
    """
    Draw the domain edges.
    """
    # draw the domain edges such that the labels are correctly positioned
    arr = []
    # draw the first edge
    line = '\\draw[thick] (n0) -- node[below] {$d_1$} (n1);'
    arr.append(line)
    # draw the second edge
    line = '\\draw[thick] (n2) -- node[below] {$d_2$} (n0);'
    arr.append(line)
    # draw the third edge
    line = '\\draw[thick] (n3) -- node[below right,inner sep=0pt] {$d_3$} (n0);'
    arr.append(line)
    return arr

def get_codomain_vertex_tikz_lines(points):
    arr = []
    for i, (x, y) in enumerate(points):
        line = '\\path (%04f, %04f) coordinate (y%s);' % (x, y, i)
        arr.append(line)
    return arr

def get_codomain_edge_tikz_lines():
    """
    Draw the four labeled gray height lines and the unlabeled black cap lines.
    """
    arr = []
    for i in range(4):
        node = 'node[black,right]'
        if i:
            label = '$y_%s$' % i
            draw = 'draw[gray]'
        else:
            label = '$y$'
            draw = 'draw[gray,dashed]'
        line = '\\%s (n%s) -- %s {%s} (y%s);' % (draw, i, node, label, i)
        arr.append(line)
    arr.extend([
        '\\draw[gray] (y0) -- (y1);',
        '\\draw[gray] (y0) -- (y2);',
        '\\draw[gray] (y0) -- (y3);'])
    return arr

def get_tikz_lines(fs):
    """
    @param fs: user input
    @return: a sequence of tikz lines
    """
    # get the harmonic extension
    w1 = 1 / fs.d1
    w2 = 1 / fs.d2
    w3 = 1 / fs.d3
    numerator = fs.y1 * w1 + fs.y2 * w2 + fs.y3 * w3
    denominator = w1 + w2 + w3
    y0 = numerator / denominator
    # get the domain point locations
    domain_points = [
            np.array([0, 0, 0], dtype=float),
            np.array([0, -fs.d1, 0], dtype=float),
            np.array([0, fs.d2, 0], dtype=float),
            np.array([-fs.d3, 0, 0], dtype=float)]
    # get the codomain point locations
    codomain_points = [np.array(p) for p in domain_points]
    codomain_points[0][2] = y0
    codomain_points[1][2] = fs.y1
    codomain_points[2][2] = fs.y2
    codomain_points[3][2] = fs.y3
    # rotate all of the points
    domain_points = [rotate_to_view(p) for p in domain_points]
    codomain_points = [rotate_to_view(p) for p in codomain_points]
    # project onto a 2D plane
    domain_2d = [p[1:] for p in domain_points]
    codomain_2d = [p[1:] for p in codomain_points]
    # get the tikz lines
    lines = []
    lines.extend([
        '% variable settings:',
        '%% d1: %04f' % fs.d1,
        '%% d2: %04f' % fs.d2,
        '%% d3: %04f' % fs.d3,
        '%% y1: %04f' % fs.y1,
        '%% y2: %04f' % fs.y2,
        '%% y3: %04f' % fs.y3])
    lines.extend(get_domain_vertex_tikz_lines(domain_2d))
    lines.extend(get_codomain_vertex_tikz_lines(codomain_2d))
    lines.extend(get_domain_edge_tikz_lines())
    lines.extend(get_codomain_edge_tikz_lines())
    return lines

def get_response_content(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: the response
    """
    tikz_body = '\n'.join(get_tikz_lines(fs))
    tikzpicture = tikz.get_picture(tikz_body, 'auto')
    return tikz.get_response(tikzpicture, fs.tikzformat)
