"""
Use tikz to draw a sample path from a continuous time Markov chain.
"""

from StringIO import StringIO

import random
import string
import numpy as np
import scipy
from scipy import linalg

import Form
import FormOut
import tikz
import latexutil
import iterutils
import mrate

g_default_matrix_string = """\
-2   2  0
49 -98 49
 0   2 -2
"""

g_default_matrix = np.array([
    [-2, 2, 0],
    [49, -98, 49],
    [0, 2, -2]])

g_default_colors = ['blue', 'red', 'yellow']


def get_form():
    """
    @return: the body of a form
    """
    # define the form objects
    form_objects = [
            Form.Matrix('R', 'rate matrix',
                g_default_matrix),
            Form.Sequence('colorseq', 'user-defined colors',
                g_default_colors),
            Form.Float('interval', 'amount of elapsed time',
                '1.0', low_exclusive=0, high_exclusive=100),
            Form.Float('tikzlength', 'tikz line length',
                '6.0', low_exclusive=1, high_exclusive=200),
            Form.TikzFormat()]
    return form_objects

def get_form_out():
    return FormOut.Tikz()

def sample_from_weights(weights):
    """
    Sample from a finite distribution defined by nonnegative weights.
    """
    probs = np.array(weights, dtype=float) / sum(weights)
    u = random.random()
    total = 0
    for i, p in enumerate(probs):
        total += p
        if u <= total:
            return i

def get_draw_line(c, x1, y1, x2, y2):
    return '\\draw[color=%s,very thick] (%0.4f,%0.4f) -- (%0.4f,%0.4f);' % (
            c, x1, y1, x2, y2)

def get_tikz_lines(fs):
    # read the user variables
    R = np.array(fs.R)
    colors = fs.colorseq
    invalid_colors = [s for s in colors if set(s) - set(string.lowercase)]
    if invalid_colors:
        raise ValueError('invalid colors: ' + str(invalid_colors))
    t_max = fs.interval
    length_total = fs.tikzlength
    # sample the initial state
    nstates = len(R)
    distn = mrate.R_to_distn(R)
    states = [sample_from_weights(distn)]
    times = [0.0]
    # sample the path
    while True:
        cur_state = states[-1]
        rate = -R[cur_state, cur_state]
        t = random.expovariate(rate)
        t_next = times[-1] + t
        if t_next >= t_max:
            break
        else:
            times.append(t_next)
            weights = np.maximum(R[cur_state], 0)
            states.append(sample_from_weights(weights))
    states.append(states[-1])
    times.append(t_max)
    # draw the path
    nsegs = len(states) - 1
    out = StringIO()
    lines = []
    for i in range(nsegs):
        state_a, state_b = states[i:i+2]
        times_a, times_b = times[i:i+2]
        xa = ((times_a / t_max) - 0.5) * length_total
        xb = ((times_b / t_max) - 0.5) * length_total
        y = 0
        line = get_draw_line(colors[state_a], xa, y, xb, y)
        lines.append(line)
    return lines


def get_response_content(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: the response
    """
    tikz_body = '\n'.join(get_tikz_lines(fs))
    tikzpicture = tikz.get_picture(tikz_body, 'auto')
    return tikz.get_response(tikzpicture, fs.tikzformat)

