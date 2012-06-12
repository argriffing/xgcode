"""
Draw two proportion sequence identity curves with different rates.

Show that there is a tradeoff between
divergence time information at small and large times.
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
                '9', low_exclusive=0, high_exclusive=20),
            Form.Float('plot_height', 'plot height in tikz units',
                '6', low_exclusive=0, high_exclusive=20),
            Form.Float('t_max', 'max time',
                '5', low_exclusive=0),
            Form.Float('slow_mu', 'slow randomization rate',
                '0.4', low_exclusive=0),
            Form.Float('fast_mu', 'fast randomization rate',
                '1', low_exclusive=0),
            Form.Float('p_width', 'proportion interval width',
                '0.06', low_exclusive=0, high_exclusive=0.75),
            Form.Float('slow_high', 'high proportion for slow process',
                '0.9', low_exclusive=0.25, high_exclusive=1),
            Form.Float('fast_high', 'high proportion for fast process',
                '0.77', low_exclusive=0.25, high_exclusive=1),
            Form.Float('slow_low', 'low proportion for slow process',
                '0.5', low_exclusive=0.25, high_exclusive=1),
            Form.Float('fast_low', 'low proportion for fast process',
                '0.3', low_exclusive=0.25, high_exclusive=1),
            #Form.Float('p_low', 'proportion lower bound',
                #'0.6', low_exclusive=0.25, high_exclusive=1.0),
            #Form.Float('p_high', 'proportion upper bound',
                #'0.8', low_exclusive=0.25, high_exclusive=1.0),
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

class Process:
    def __init__(self,
            plot_width, plot_height, timescale, p_width,
            mu, p_low, p_high):
        self.plot_width = plot_width
        self.plot_height = plot_height
        self.timescale = timescale
        self.p_width = p_width
        self.mu = mu
        self.p_low = p_low
        self.p_high = p_high
        # validate
        if p_high <= p_low:
            raise ValueError(
                    'interval lower bound should be below upper bound')
        # aux members
        self.f = MyCurve(self.mu)
    def _get_knot_times(self):
        return sorted((
                0.0, self.timescale,
                self.f.inv(self.p_low - 0.5*self.p_width),
                self.f.inv(self.p_low + 0.5*self.p_width),
                self.f.inv(self.p_high - 0.5*self.p_width),
                self.f.inv(self.p_high + 0.5*self.p_width)))
    def draw_curve(self):
        scale = np.array((self.plot_width / self.timescale, self.plot_height))
        times = self._get_knot_times()
        bchunks = []
        for a, b in iterutils.pairwise(times):
            pta = np.array((a, self.f(a)))
            ptb = np.array((b, self.f(b)))
            dta = np.array((1, self.f.deriv(a)))
            dtb = np.array((1, self.f.deriv(b)))
            bchunk = bezier.create_bchunk_hermite(
                    a, b,
                    pta * scale, ptb * scale,
                    dta * scale, dtb * scale)
            bchunks.append(bchunk)
        return r'\draw ' + get_tikz_bezier(bchunks)
    def _draw_beam(self, p_mid, color):
        scale = np.array((self.plot_width / self.timescale, self.plot_height))
        xproj = np.array((1, 0))
        yproj = np.array((0, 1))
        out = StringIO()
        print >> out, r'\path[fill=%s,fill opacity=0.5]' % color
        p_upper = p_mid + 0.5*self.p_width
        p_lower = p_mid - 0.5*self.p_width
        t_upper = self.f.inv(p_lower)
        t_lower = self.f.inv(p_upper)
        pta = np.array((t_lower, p_upper))
        ptb = np.array((t_upper, p_lower))
        dta = np.array((1, self.f.deriv(t_lower)))
        dtb = np.array((1, self.f.deriv(t_upper)))
        print >> out, tikz.point_to_tikz(pta*scale*yproj) + ' --'
        bchunk = bezier.create_bchunk_hermite(
                t_lower, t_upper,
                pta * scale, ptb * scale,
                dta * scale, dtb * scale)
        pts = tuple(tikz.point_to_tikz(p) for p in bchunk.get_points())
        print >> out, '%s .. controls %s and %s .. %s --' % pts
        print >> out, tikz.point_to_tikz(ptb*scale*xproj) + ' --'
        print >> out, tikz.point_to_tikz(pta*scale*xproj) + ' --'
        ptc = np.array((t_lower, p_lower))
        print >> out, tikz.point_to_tikz(ptc*scale) + ' --'
        print >> out, tikz.point_to_tikz(ptb*scale*yproj) + ' -- cycle;'
        return out.getvalue().rstrip()
    def draw_high_beam(self, color):
        return self._draw_beam(self.p_high, color)
    def draw_low_beam(self, color):
        return self._draw_beam(self.p_low, color)


def get_tikz_body(fs):
    out = StringIO()
    # init the processes from user data
    fast_process = Process(
            fs.plot_width, fs.plot_height, fs.t_max, fs.p_width,
            fs.fast_mu, fs.fast_low, fs.fast_high)
    slow_process = Process(
            fs.plot_width, fs.plot_height, fs.t_max, fs.p_width,
            fs.slow_mu, fs.slow_low, fs.slow_high)
    # predefined variables
    origin = (0, 0)
    # define user variables
    plot_width = fs.plot_width
    plot_height = fs.plot_height
    timescale = fs.t_max
    plotscale = np.array((plot_width, plot_height))
    # draw the beams
    print >> out, slow_process.draw_high_beam('blue!60')
    print >> out, slow_process.draw_low_beam('blue!60')
    print >> out, fast_process.draw_high_beam('red!60')
    print >> out, fast_process.draw_low_beam('red!60')
    # draw the boundary of the plot
    print >> out, r'\draw[color=gray] ' + get_segment(
            origin, (plot_width, 0))
    print >> out, r'\draw[color=gray] ' + get_segment(
            origin, (0, plot_height))
    print >> out, r'\draw[color=gray] ' + get_segment(
            (0,plot_height), (plot_width, plot_height))
    print >> out, r'\draw[dotted,color=gray] ' + get_segment(
            (0,0.25*plot_height), (plot_width, 0.25*plot_height))
    # draw the bezier curves hitting the right knots
    for p in (slow_process, fast_process):
        print >> out, p.draw_curve()
    # draw filled black dots at some intersections
    dot_points = [
            origin,
            (0, plot_height),
            (0, 0.25 * plot_height),
            ]
    for dot_point in dot_points:
        print >> out, r'\fill[color=black,inner sep=0pt]',
        print >> out, tikz.point_to_tikz(dot_point),
        print >> out, 'circle (1pt);'
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
    return tikz.get_response(tikzpicture, fs.tikzformat)

