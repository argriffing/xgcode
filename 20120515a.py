"""
Draw proportion sequence identity curve as a function of time.

The specific model should not really matter
because this is a cartoon.
Jukes-Cantor could be a good model
because a proportion identity asymptote at 0.25 might look better than others.
Try using a bezier to draw the curve.
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

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.Float('plot_width', 'plot width in tikz units',
                '4', low_exclusive=0, high_exclusive=20),
            Form.Float('plot_height', 'plot height in tikz units',
                '6', low_exclusive=0, high_exclusive=20),
            Form.Float('t_max', 'max time',
                '5', low_exclusive=0),
            Form.Float('p_low', 'proportion lower bound',
                '0.6', low_exclusive=0.25, high_exclusive=1.0),
            Form.Float('p_high', 'proportion upper bound',
                '0.8', low_exclusive=0.25, high_exclusive=1.0),
            #Form.Float('t_a', 'lower bound of interval of interest',
                #'1.0', low_exclusive=0),
            #Form.Float('t_b', 'upper bound of interval of interest',
                #'1.5', low_exclusive=0),
            Form.TikzFormat()]
    return form_objects

def get_form_out():
    return FormOut.Tikz()

class MyCurve:
    def __init__(self, mu):
        """
        This is P(X(0) == X(t)) for 4-state Jukes-Cantor.
        @param mu: randomization rate
        """
        self.mu = mu
        # define the logical entropy of the stationary distribution
        self.h = 0.75
    def deriv(self, t):
        return -self.h * self.mu * math.exp(-self.mu * t)
    def inv(self, p):
        return -math.log((p + self.h - 1) / self.h) / self.mu
    def __call__(self, t):
        return self.h * math.exp(-self.mu * t) + (1 - self.h)

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
    # predefined variables
    mu = 1.0
    origin = (0, 0)
    f = MyCurve(mu)
    # define user variables
    plot_width = fs.plot_width
    plot_height = fs.plot_height
    timescale = fs.t_max
    ta = f.inv(fs.p_high) / timescale
    tb = f.inv(fs.p_low) / timescale
    # validate
    if tb <= ta:
        raise ValueError('interval lower bound should be below upper bound')
    plotscale = np.array((plot_width, plot_height))
    # draw the boundary of the plot
    print >> out, r'\draw[color=gray] ' + get_segment(
            origin, (plot_width, 0))
    print >> out, r'\draw[color=gray] ' + get_segment(
            origin, (0, plot_height))
    print >> out, r'\draw[color=gray] ' + get_segment(
            (0,plot_height), (plot_width, plot_height))
    print >> out, r'\draw[dotted,color=gray] ' + get_segment(
            (0,0.25*plot_height), (plot_width, 0.25*plot_height))
    # define times of interest
    t0 = 0
    tx = (tb + 1) / 2
    t1 = 1
    # draw the bezier curve hitting the right knots
    scale = np.array((plot_width / timescale, plot_height))
    times = (t0, ta, tb, tx, t1)
    bchunks = []
    for a, b in iterutils.pairwise(times):
        a = timescale * a
        b = timescale * b
        pta = np.array((a, f(a)))
        ptb = np.array((b, f(b)))
        dta = np.array((1, f.deriv(a)))
        dtb = np.array((1, f.deriv(b)))
        bchunk = bezier.create_bchunk_hermite(
                a, b,
                pta * scale, ptb * scale,
                dta * scale, dtb * scale)
        bchunks.append(bchunk)
    print >> out, r'\draw[color=gray] ' + get_tikz_bezier(bchunks)
    # redraw a piece of the curve
    a = timescale * ta
    b = timescale * tb
    pta = np.array((a, f(a)))
    ptb = np.array((b, f(b)))
    dta = np.array((1, f.deriv(a)))
    dtb = np.array((1, f.deriv(b)))
    bchunk = bezier.create_bchunk_hermite(
            a, b,
            pta * scale, ptb * scale,
            dta * scale, dtb * scale)
    pts = tuple(tikz.point_to_tikz(p) for p in bchunk.get_points())
    print >> out, r'\draw[color=black] %s .. controls %s and %s .. %s;' % pts
    """
    print >> out, r'\draw[color=black] ' + get_segment(
            pta * scale, ptb * scale)
    """
    # draw the projections of the secant onto the axes
    xproj = np.array((1, 0))
    yproj = np.array((0, 1))
    print >> out, r'\draw[color=black] ' + get_segment(
            pta * scale * xproj, ptb * scale * xproj)
    print >> out, r'\draw[color=black] ' + get_segment(
            pta * scale * yproj, ptb * scale * yproj)
    print >> out, r'\draw[dotted,color=gray] ' + get_segment(
            pta * scale, pta * scale * xproj)
    print >> out, r'\draw[dotted,color=gray] ' + get_segment(
            ptb * scale, ptb * scale * xproj)
    print >> out, r'\draw[dotted,color=gray] ' + get_segment(
            pta * scale, pta * scale * yproj)
    print >> out, r'\draw[dotted,color=gray] ' + get_segment(
            ptb * scale, ptb * scale * yproj)
    # draw filled black dots at some intersections
    dot_points = [
            origin,
            (0, plot_height),
            (0, 0.25 * plot_height),
            pta * scale,
            ptb * scale,
            pta * scale * xproj,
            pta * scale * yproj,
            ptb * scale * xproj,
            ptb * scale * yproj,
            ]
    for dot_point in dot_points:
        print >> out, r'\fill[color=black,inner sep=0pt]',
        print >> out, tikz.point_to_tikz(dot_point),
        print >> out, 'circle (1pt);'
    # draw braces
    brace_terms = [
        r'\draw[decorate,decoration={brace},yshift=-2pt] ',
        get_seg(ptb * scale * xproj, pta * scale * xproj),
        r'node [black,midway,yshift=-2pt]',
        #r'{$\Delta t_{\text{divergence}}$};']
        r'{$\Delta t$};']
    print >> out, ' '.join(brace_terms)
    brace_terms = [
        r'\draw[decorate,decoration={brace},xshift=-2pt] ',
        get_seg(ptb * scale * yproj, pta * scale * yproj),
        r'node [black,midway,xshift=-2pt]',
        #r'{$\Delta P_{\text{identity}}$};']
        r'{$\Delta p$};']
    print >> out, ' '.join(brace_terms)
    # draw some text annotations
    pt_txt_pairs = [
            ((0, 0), '0'),
            ((0, 0.25 * plot_height), r'$\frac{1}{4}$'),
            ((0, 1.0 * plot_height), '1')]
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
    tikzlibraries = ['decorations.pathreplacing']
    return tikz.get_response(
            tikzpicture, fs.tikzformat, tikzlibraries=tikzlibraries)

