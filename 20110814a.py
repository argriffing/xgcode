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

def get_tikz_pane(sample, width=6, height=6):
    abstol = 1e-6
    vgap = 0.5
    cut_radius = 0.1
    shape = sample.get_shape()
    try:
        time_lists = shape.get_orthoplanar_intersection_times()
        t_initial = shape.t_initial
        t_final = shape.t_final
    except AttributeError, e:
        return '\\node (0, 0) {sign cut not available};'
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

def get_latex_text(tikz_text):
    """
    TikZ boilerplate code.
    """
    preamble = '\\usepackage{color}'
    arr = [tikz.define_color(*p) for p in color.wolfram_name_color_pairs]
    arr.append(tikz_text)
    document_body = '\n'.join(arr)
    return tikz.get_latex_text(preamble, document_body)

def get_tikz_text(tikz_body):
    """
    TikZ boilerplate code.
    """
    return '\n'.join([
            '\\begin{tikzpicture}[auto]',
            tikz_body,
            '\\end{tikzpicture}'])

def get_response_content(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: the response
    """
    # get the texts
    tikz_lines = get_tikz_lines()
    tikz_text = get_tikz_text('\n'.join(tikz_lines))
    latex_text = get_latex_text(tikz_text)
    # decide the output format
    if fs.tikz:
        return tikz_text
    elif fs.tex:
        return latex_text
    elif fs.pdf:
        return tikz.get_pdf_contents(latex_text)
    elif fs.png:
        return tikz.get_png_contents(latex_text)


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

