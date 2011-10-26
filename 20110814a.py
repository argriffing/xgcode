"""
Draw nine hardcoded sign cut figures.

The command line version will create tikz files.
"""

from StringIO import StringIO
import argparse
import math
import os
import sys

import numpy as np

import Form
import FormOut
import tikz
import interlace
import interlacesample
import pcurve
import sympyutils
import color
import typeutils
import iterutils


g_colors = ['black'] + color.wolfram[:3]

def get_form():
    """
    @return: a list of form objects
    """
    # define the form objects
    form_objects = [
            Form.TikzFormat(),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.Tikz()

def get_tree_tikz_pane(
        shape, width, height, time_lists,
        t_initial, t_final,
        vgap, cut_radius):
    pass

def get_linear_tikz_pane(
        shape, width, height, time_lists,
        t_initial, t_final,
        vgap, cut_radius):
    abstol = 1e-6
    duration = float(t_final - t_initial)
    arr = []
    for i in range(4):
        c = g_colors[i]
        xa = t_initial * (width / duration)
        xb = t_final * (width / duration)
        # draw the thin line of the correct color
        line = '\\draw[%s] %s -- %s;' % (
                c,
                tikz.point_to_tikz((xa, -i*vgap)),
                tikz.point_to_tikz((xb, -i*vgap)))
        arr.append(line)
        # draw the thick segments of the correct color
        if i:
            augmented_times = [t_initial] + time_lists[i-1] + [t_final]
            for ta, tb in iterutils.pairwise(augmented_times):
                t = (ta + tb) / 2.0
                xa = ta * (width / duration)
                xb = tb * (width / duration)
                value = shape.fps[i-1](t)
                if value > 0:
                    line = '\\draw[very thick,%s] %s -- %s;' % (
                            c,
                            tikz.point_to_tikz((xa, -i*vgap)),
                            tikz.point_to_tikz((xb, -i*vgap)))
                    arr.append(line)
        # draw the cuts in black ink
        if i < 3:
            times = time_lists[i]
            for t in times:
                x = t * (width / duration)
                line = '\\draw %s -- %s;' % (
                        tikz.point_to_tikz((x, cut_radius-i*vgap)),
                        tikz.point_to_tikz((x, -cut_radius-i*vgap)))
                arr.append(line)
    return '\n'.join(arr)

def get_tree_tikz_pane(sample, width, height, vgap, cut_radius):
    # Use the target width to scale the layout.
    # TODO rotate the layout for greatest width to height ratio.
    unscaled_layout_points = sample.v_to_layout_point.values()
    pmin = np.min(unscaled_layout_points, axis=0)
    pmax = np.max(unscaled_layout_points, axis=0)
    sf_width = width / (pmax[0] - pmin[0])
    sf_height = height / ((pmax[1] - pmin[1]) * 4 + vgap * 3)
    sf = min(sf_width, sf_height)
    v_to_layout_point = dict(
            (v, p*sf) for v, p in sample.v_to_layout_point.items())
    # Compute the amount (in tikz units) to skip per tree.
    vskip = (pmax[1] - pmin[1])*sf + vgap
    # Draw the tikz.
    arr = []
    vskip_accum = 0.0
    for i in range(4):
        # define the current offset for drawing
        offset = np.array([0, -vskip_accum])
        # define the color of the tree
        c = g_colors[i]
        # draw the tree using thin line segments
        for va, vb in sample.T:
            pa, pb = v_to_layout_point[va], v_to_layout_point[vb]
            # draw the thin line segment of the correct color
            line = '\\draw[%s] %s -- %s;' % (
                    c,
                    tikz.point_to_tikz(pa + offset),
                    tikz.point_to_tikz(pb + offset))
            arr.append(line)
        # draw the thick segments of the correct color
        if i:
            axis = i-1
            for va, vb in sample.T:
                vala = sample.v_to_point[va][axis]
                valb = sample.v_to_point[vb][axis]
                pa, pb = v_to_layout_point[va], v_to_layout_point[vb]
                if vala > 0 and valb > 0:
                    # If both endpoints are positive
                    # then redraw the whole segment using a very thick line.
                    line = '\\draw[very thick,%s] %s -- %s;' % (
                            c,
                            tikz.point_to_tikz(pa + offset),
                            tikz.point_to_tikz(pb + offset))
                    arr.append(line)
                elif vala * valb < 0:
                    # If the endpoints have opposite sign
                    # then redraw only part of the line.
                    t = vala / (vala - valb)
                    p = (1 - t) * pa + t * pb
                    if vala > 0:
                        p_begin = pa
                        p_end = p
                    else:
                        b_begin = p
                        p_end = pb
                    line = '\\draw[very thick,%s] %s -- %s;' % (
                            c,
                            tikz.point_to_tikz(p_begin + offset),
                            tikz.point_to_tikz(p_end + offset))
                    arr.append(line)
        #TODO draw tick marks
        if i < 3:
            axis = i
            for va, vb in sample.T:
                vala = sample.v_to_point[va][axis]
                valb = sample.v_to_point[vb][axis]
                pa, pb = v_to_layout_point[va], v_to_layout_point[vb]
                if vala * valb < 0:
                    # If the endpoints have opposite sign
                    # then draw a tick mark.
                    t = vala / (vala - valb)
                    p = (1 - t) * pa + t * pb
                    theta = math.atan2((pb-pa)[1], (pb-pa)[0])
                    phi = theta + math.pi/2
                    cuta = p + offset
                    cutb = p + offset
                    cuta[0] += cut_radius * math.cos(phi)
                    cuta[1] += cut_radius * math.sin(phi)
                    cutb[0] -= cut_radius * math.cos(phi)
                    cutb[1] -= cut_radius * math.sin(phi)
                    line = '\\draw %s -- %s;' % (
                            tikz.point_to_tikz(cuta),
                            tikz.point_to_tikz(cutb))
                    arr.append(line)
        # skip some space
        vskip_accum += vskip
    # return the tikz text
    return '\n'.join(arr)

def get_tikz_pane(sample, width=6, height=6):
    vgap = 0.5
    cut_radius = 0.1
    time_lists = None
    try:
        shape = sample.get_shape()
        time_lists = shape.get_orthoplanar_intersection_times()
        t_initial = shape.t_initial
        t_final = shape.t_final
    except AttributeError as e:
        pass
    if time_lists:
        return get_linear_tikz_pane(
                shape, width, height, time_lists,
                t_initial, t_final,
                vgap, cut_radius)
    else:
        return get_tree_tikz_pane(sample, width, height, vgap, cut_radius)

def get_tikz_lines():
    """
    @param fs: user input
    @return: a sequence of tikz lines
    """
    arr = []
    # draw the matrix
    samples = interlacesample.get_samples()
    arr.extend([
        '\\matrix{',
        get_tikz_pane(samples[0]),
        '&',
        get_tikz_pane(samples[1]),
        '&',
        get_tikz_pane(samples[2]),
        '\\\\',
        get_tikz_pane(samples[3]),
        '&',
        get_tikz_pane(samples[4]),
        '&',
        get_tikz_pane(samples[5]),
        '\\\\',
        get_tikz_pane(samples[6]),
        '&',
        get_tikz_pane(samples[7]),
        '&',
        get_tikz_pane(samples[8]),
        '\\\\};'])
    return arr

def get_response_content(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: the response
    """
    tikz_body = '\n'.join(get_tikz_lines())
    tikzpicture = tikz.get_picture(tikz_body, 'auto')
    return tikz.get_response(
            tikzpicture, fs.tikzformat,
            tikz.get_w_color_package_set(), tikz.get_w_color_preamble())

def main(args):
    for i, sample in enumerate(interlacesample.get_samples()):
        filename = os.path.join(args.outdir, 'sample-%04d.tikz' % i)
        with open(filename, 'w') as fout:
            print 'writing', filename
            arr = []
            # add the remark about the invocation of the generating script
            arr.append('% ' + ' '.join(sys.argv))
            # add the commands to define custom colors
            for name, rgb in color.wolfram_name_color_pairs:
                arr.append(tikz.define_color(name, rgb))
            # add the tikz drawing functions
            arr.append(get_tikz_pane(sample))
            # write the file
            print >> fout, '\n'.join(arr)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir',
            default='', help='output directory')
    main(parser.parse_args())

