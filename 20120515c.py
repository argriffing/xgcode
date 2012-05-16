"""
Draw two derivatives of sequence identity curves with different rates.

Show that there is a tradeoff between
divergence time information at small and large times.
The Fisher information approaches infinity at time zero,
and so does the negative of the derivative of the mutual information.
"""

from StringIO import StringIO
import string
import math

import numpy as np

import Form
import FormOut
import tikz
import latexutil
import iterutils
import mrate
import bezier
import JC69

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.Float('plot_width', 'plot width in tikz units',
                '4', low_exclusive=0, high_exclusive=20),
            Form.Float('plot_height', 'plot height in tikz units',
                '4', low_exclusive=0, high_exclusive=20),
            Form.Float('t_max', 'max time',
                '5', low_exclusive=0),
            Form.Float('slow_mu', 'slow randomization rate',
                '0.4', low_exclusive=0),
            Form.Float('fast_mu', 'fast randomization rate',
                '1', low_exclusive=0),
            Form.RadioGroup('info_type', 'information type', [
                Form.RadioItem('info_identity_slope',
                    'abs slope of expected identity', True),
                Form.RadioItem('info_fi', 'Fisher information'),
                Form.RadioItem('info_mi', 'mutual information')]),
            Form.TikzFormat()]
    return form_objects

def get_form_out():
    return FormOut.Tikz()

def get_tikz_bezier(bchunks):
    """
    @param bchunks: a sequence of 2d bezier chunks
    @return: multiline bezier text
    """
    lines = []
    # draw everything except for the last point of the last chunk
    for b in bchunks:
        pts = [tikz.point_to_tikz(p) for p in b.get_points()[:-1]]
        lines.append('%s .. controls %s and %s ..' % tuple(pts))
    # draw the last point of the last chunk
    lines.append('%s;' % tikz.point_to_tikz(bchunks[-1].p3))
    return '\n'.join(lines)

def get_seg(pta, ptb):
    return '%s -- %s' % (tikz.point_to_tikz(pta), tikz.point_to_tikz(ptb))

def get_segment(pta, ptb):
    return get_seg(pta, ptb) + ';'

def get_tikz_body(fs):
    out = StringIO()
    # define user variables
    plot_width = fs.plot_width
    plot_height = fs.plot_height
    timescale = fs.t_max
    fast_mu = fs.fast_mu
    slow_mu = fs.slow_mu
    if fs.info_identity_slope:
        f_fast = JC69.IdentitySlopeInformation(fast_mu)
        f_slow = JC69.IdentitySlopeInformation(slow_mu)
    elif fs.info_mi:
        f_fast = JC69.MutualInformation(fast_mu)
        f_slow = JC69.MutualInformation(slow_mu)
    elif fs.info_fi:
        #f_fast = JC69.FisherInformationTheano(fast_mu)
        #f_slow = JC69.FisherInformationTheano(slow_mu)
        f_fast = JC69.FisherInformation(fast_mu)
        f_slow = JC69.FisherInformation(slow_mu)
    # Define some times for evaluation of the curve.
    times = [timescale*2**-i for i in range(10)]
    if fs.info_identity_slope:
        # Compute the intersection time.
        t_x = math.log(fast_mu / slow_mu) / (fast_mu - slow_mu)
        times.extend([t_x / 2, t_x, (t_x + timescale)/2])
    # define some more intermediate values
    ymax = max(f_fast(min(times)), f_slow(min(times))) * 1.2
    plotscale = np.array((plot_width / timescale, plot_height / ymax))
    origin = (0, 0)
    # draw the boundary of the plot
    print >> out, r'\draw[color=gray] %s %s {%s} %s;' % (
            tikz.point_to_tikz(origin),
            'edge node[color=black,below]',
            '$t$',
            tikz.point_to_tikz((plot_width, 0)))
    print >> out, r'\draw[color=gray] ' + get_segment(
            origin, (0, plot_height))
    # draw the bezier curves hitting the right knots
    for f in (f_slow, f_fast):
        bchunks = []
        for a, b in iterutils.pairwise(sorted(times)):
            pta = np.array((a, f(a)))
            ptb = np.array((b, f(b)))
            dta = np.array((1, f.deriv(a)))
            dtb = np.array((1, f.deriv(b)))
            bchunk = bezier.create_bchunk_hermite(
                    a, b,
                    pta * plotscale, ptb * plotscale,
                    dta * plotscale, dtb * plotscale)
            bchunks.append(bchunk)
        print >> out, r'\draw[color=gray] ' + get_tikz_bezier(bchunks)
    # draw filled black dots at some intersections
    dot_points = [origin]
    if not fs.info_fi:
        dot_points.append((0, f_fast(0)))
        dot_points.append((0, f_slow(0)))
    if fs.info_identity_slope:
        dot_points.append((t_x, f_slow(t_x)))
    for p in dot_points:
        print >> out, r'\fill[color=black,inner sep=0pt]',
        print >> out, tikz.point_to_tikz(np.array(p) * plotscale),
        print >> out, 'circle (1pt);'
    # draw some text annotations
    pt_txt_pairs = [
            ((0, 0), '0'),
            ]
    for i, (pt, txt) in enumerate(pt_txt_pairs):
        print >> out, r'\node[anchor=east] (%s) at %s {%s};' % (
                'ylabel%d' % i,
                tikz.point_to_tikz(pt),
                txt)
    #
    return out.getvalue().rstrip()


def get_response_content(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: the response
    """
    tikz_body = get_tikz_body(fs)
    tikzpicture = tikz.get_picture(tikz_body, 'auto')
    return tikz.get_response(tikzpicture, fs.tikzformat)

