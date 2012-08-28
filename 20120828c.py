"""
Reproduce yet another figure from the paper by Christina Chen et al.
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
    low = -20
    high = 20
    s_values = np.linspace(low, high, (high-low)*10 + 1)
    # define some dominance coefficients
    h_values = [-3, 0, 1, 2]
    # get the values for each h
    arr = []
    for h in h_values:
        v = [kimura.get_fixation_probability_chen(p, s, h) for s in s_values]
        arr.append(v)
    #
    # define the r script
    out = StringIO()
    print >> out, 's.values <- c', str(tuple(s_values))
    print >> out, 'ha <- c', str(tuple(arr[0]))
    print >> out, 'hb <- c', str(tuple(arr[1]))
    print >> out, 'hc <- c', str(tuple(arr[2]))
    print >> out, 'hd <- c', str(tuple(arr[3]))
    print >> out, mk_call_str('plot', 's.values', 'ha',
            type='"l"',
            xlab='"selection coefficient (s)"',
            ylab='"probability of eventual fixation"',
            main='"fixation probabilities for various h ; p0 = 0.1"',
            ylim='c(0,1)',
            )
    print >> out, mk_call_str('lines', 's.values', 'hb', col='"red"')
    print >> out, mk_call_str('lines', 's.values', 'hc', col='"green"')
    print >> out, mk_call_str('lines', 's.values', 'hd', col='"blue"')
    script = out.getvalue().rstrip()
    # create the R plot image
    device_name = Form.g_imageformat_to_r_function[fs.imageformat]
    retcode, r_out, r_err, image_data = RUtil.run_plotter_no_table(
            script, device_name)
    if retcode:
        raise RUtil.RError(r_err)
    return image_data

