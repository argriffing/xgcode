"""
Check whether dominance and selection are confounded near h=1/2.
"""

from StringIO import StringIO
import math

import numpy as np
from scipy import special

import Form
import FormOut
import RUtil
from RUtil import mk_call_str
import kimura

def get_form():
    """
    @return: the body of a form
    """
    return [
            Form.ImageFormat(),
            ]

def get_form_out():
    return FormOut.Image('plot')

def get_response_content(fs):
    # define the initial frequency
    p = 0.1
    # define some selection coefficients to plot
    h_low = 0.2
    h_high = 0.8
    h_values = np.linspace(h_low, h_high, 6*10 + 1)
    # define some dominance coefficients
    hs_values = [-0.2, -0.05, 0, 0.05, 0.2]
    colors = ['blue', 'green', 'black', 'orange', 'red']
    # get the values for each h
    arr = []
    for hs in hs_values:
        v = [kimura.get_fixation_probability_chen(p, hs/h, h) for h in h_values]
        arr.append(v)
    #
    # define the r script
    out = StringIO()
    print >> out, 'h.values <- c', str(tuple(h_values))
    print >> out, 'ha <- c', str(tuple(arr[0]))
    print >> out, 'hb <- c', str(tuple(arr[1]))
    print >> out, 'hc <- c', str(tuple(arr[2]))
    print >> out, 'hd <- c', str(tuple(arr[3]))
    print >> out, 'he <- c', str(tuple(arr[4]))
    print >> out, mk_call_str('plot', 'h.values', 'ha',
            type='"l"',
            xlab='"h"',
            ylab='"probability of eventual fixation"',
            main=(
                '"fixation probabilities for various h*s ; '
                's = 2*N*sigma ; p0 = %s"' % p),
            ylim='c(0, 0.2)',
            col='"%s"' % colors[0],
            )
    print >> out, mk_call_str('lines', 'h.values', 'hb', col='"%s"' % colors[1])
    print >> out, mk_call_str('lines', 'h.values', 'hc', col='"%s"' % colors[2])
    print >> out, mk_call_str('lines', 'h.values', 'hd', col='"%s"' % colors[3])
    print >> out, mk_call_str('lines', 'h.values', 'he', col='"%s"' % colors[4])
    script = out.getvalue().rstrip()
    # create the R plot image
    device_name = Form.g_imageformat_to_r_function[fs.imageformat]
    retcode, r_out, r_err, image_data = RUtil.run_plotter_no_table(
            script, device_name)
    if retcode:
        raise RUtil.RError(r_err)
    return image_data

