"""
Draw Bezier curve intersection. [BIT ROTTED]
"""

import math

import numpy as np
import sympy

import Form
import FormOut
import tikz
import color
import pcurve
import bezier

def get_form():
    """
    @return: a list of form objects
    """
    # define the form objects
    form_objects = [
            Form.Float('min_gridsize', 'intersection resolution',
                0.1, low_exclusive=0),
            Form.TikzFormat()]
    return form_objects

def get_form_out():
    return FormOut.Tikz()

def get_tikz_lines(fs):
    """
    @param fs: user input
    @return: a sequence of tikz lines
    """
    # define a straightforward bezier curve
    p0 = np.array([-1.0, -1.0])
    p1 = np.array([-1.0, 1.0])
    p2 = np.array([1.0, 1.0])
    p3 = np.array([1.0, -1.0])
    # define a messier bezier curve
    q0 = np.array([-2.0, 0.0])
    q1 = np.array([1.0, 1.0])
    q2 = np.array([-1.0, -1.0])
    q3 = np.array([2.0, 0.0])
    # plot the bezier curves using tikz
    arr = []
    points = tuple(tikz.point_to_tikz(p) for p in (p0, p1, p2, p3))
    arr.append('\\draw %s .. controls %s and %s .. %s;' % points)
    points = tuple(tikz.point_to_tikz(p) for p in (q0, q1, q2, q3))
    arr.append('\\draw %s .. controls %s and %s .. %s;' % points)
    #
    a = bezier.BezierChunk()
    a.p0 = p0
    a.p1 = p1
    a.p2 = p2
    a.p3 = p3
    a.start_time = 0.0
    a.stop_time = 1.0
    a.parent_ref = 10
    #
    b = bezier.BezierChunk()
    b.p0 = q0
    b.p1 = q1
    b.p2 = q2
    b.p3 = q3
    b.start_time = 0.0
    b.stop_time = 1.0
    b.parent_ref = 11
    # find the intersections
    beziers = pcurve.find_bezier_intersections([a, b], fs.min_gridsize)
    for b in beziers:
        points = tuple(tikz.point_to_tikz(p) for p in (b.p0, b.p1, b.p2, b.p3))
        arr.append('\\draw[red] %s .. controls %s and %s .. %s;' % points)
    # return the lines
    return arr

def get_response_content(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: the response
    """
    tikz_body = '\n'.join(get_tikz_lines(fs))
    tikzpicture = tikz.get_picture(tikz_body, 'auto')
    return tikz.get_response(
            tikzpicture, fs.tikzformat,
            tikz.get_w_color_package_set(), tikz.get_w_color_preamble())
