"""
Check whether dominance and selection are confounded near h=1/2, even more.
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
            #Form.FloatInterval(
                #'x_min', 'x_max', 'log probability ratio interval', '-1', '2'),
            Form.Float('s_scale', 'scale the selection interval',
                1.0, low_inclusive=0.01, high_inclusive=10.0),
            Form.RadioGroup('linetype', 'h values', [
                Form.RadioItem('h_dense', '0.2 < h < 0.8', True),
                Form.RadioItem('h_sparse', 'h in {0, 0.5, 1}'),
                ]),
            Form.RadioGroup('xoptions', 'x axis', [
                Form.RadioItem(
                    'x_is_pfix', 'fixation probability', True),
                Form.RadioItem(
                    'x_is_log_pfix', 'log of fixation probability ratio'),
                ]),
            Form.RadioGroup('yoptions', 'y axis', [
                Form.RadioItem('infer_s', 's', True),
                Form.RadioItem('infer_hs', 'hs'),
                ]),
            Form.ImageFormat(),
            ]

def get_form_out():
    return FormOut.Image('plot')

"""
def get_presets():
    return [
            Form.Preset('small selection', {'s_scale' : '0.10'}),
            Form.Preset('large selection', {'s_scale' : '10.0'}),
            ]
"""

def rev(x):
    return tuple(reversed(list(x)))

def get_response_content(fs):
    # define the initial frequency
    p = 0.01
    density = 2
    # dominance
    if fs.h_dense:
        h_low = 0.2
        h_high = 0.8
        h_values = np.linspace(h_low, h_high, 6*density + 1)
    elif fs.h_sparse:
        h_low = 0
        h_high = 1
        h_values = np.linspace(h_low, h_high, 2 + 1)
    # selection
    s_low = -1.0
    s_high = 1.0
    s_values = fs.s_scale * np.linspace(s_low, s_high, 20*density + 1)
    # hs product
    hs_values = np.outer(h_values, s_values)
    # compute all values in the grid
    f_values = np.zeros((len(h_values), len(s_values)))
    for i, h in enumerate(h_values):
        for j, s in enumerate(s_values):
            f_values[i, j] = kimura.get_fixation_probability_chen(p, s, h)
    if fs.x_is_log_pfix:
        f_values = np.log(f_values / p)
    # define the r script
    out = StringIO()
    if fs.infer_s:
        ylab_string = '"selection coefficient (s)"'
        ylim_low = np.min(s_values)
        ylim_high = np.max(s_values)
    elif fs.infer_hs:
        ylab_string = '"product of dominance and selection (hs)"'
        ylim_low = np.min(hs_values)
        ylim_high = np.max(hs_values)
    if fs.x_is_pfix:
        xlab_string = '"probability of eventual fixation"'
    elif fs.x_is_log_pfix:
        xlab_string = '"log of ratios of probability of eventual fixation"'
    xlim_low = np.min(f_values)
    xlim_high = np.max(f_values)
    colors = ['orange', 'green', 'blue']
    if fs.h_dense:
        main_string = (
                '"s = 2*N*sigma ; '
                '%s < h < %s ; '
                'p0 = %s"' % (h_low, h_high, p))
    elif fs.h_sparse:
        main_string = (
                '"s = 2*N*sigma ; '
                'p0 = %s"' % p)
    print >> out, mk_call_str(
            'plot',
            'c(0,0)',
            'c(0,0)',
            type='"n"',
            xlab=xlab_string,
            ylab=ylab_string,
            main=main_string,
            xlim='c(%s, %s)' % (xlim_low, xlim_high),
            ylim='c(%s, %s)' % (ylim_low, ylim_high),
            )
    if fs.h_dense:
        for h, s_row in zip(h_values, f_values):
            if fs.infer_s:
                y_vect = s_values
            elif fs.infer_hs:
                y_vect = h*s_values
            print >> out, mk_call_str(
                    'lines',
                    'c' + str(tuple(s_row)),
                    'c' + str(tuple(y_vect)),
                    )
        for s, h_col in zip(s_values, f_values.T):
            if fs.infer_s:
                y_vect = s*np.ones_like(h_col)
            elif fs.infer_hs:
                y_vect = s*h_values
            print >> out, mk_call_str(
                    'lines',
                    'c' + str(tuple(h_col)),
                    'c' + str(tuple(y_vect)),
                    )
    elif fs.h_sparse:
        for i, (h, s_row) in enumerate(zip(h_values, f_values)):
            if fs.infer_s:
                y_vect = s_values
            elif fs.infer_hs:
                y_vect = h*s_values
            print >> out, mk_call_str(
                    'lines',
                    'c' + str(tuple(s_row)),
                    'c' + str(tuple(y_vect)),
                    col='"%s"' % colors[i],
                    )
    if fs.h_sparse:
        print >> out, mk_call_str(
                'legend',
                '"topleft"',
                'c' + str(tuple('h = %s' % x for x in h_values)),
                lty='c' + str(tuple([1]*len(h_values))),
                lwd='c' + str(tuple([2.5]*len(h_values))),
                col='c' + str(tuple(colors)),
                )
    script = out.getvalue().rstrip()
    # create the R plot image
    device_name = Form.g_imageformat_to_r_function[fs.imageformat]
    retcode, r_out, r_err, image_data = RUtil.run_plotter_no_table(
            script, device_name)
    if retcode:
        raise RUtil.RError(r_err)
    return image_data

