"""
Plot expected log likelihood ratios over time.
"""

from StringIO import StringIO
import argparse
import math

import numpy as np
import scipy
from scipy import linalg

import Form
import FormOut
from MatrixUtil import ndot
import mrate
import RUtil

UNIFORM = 'uniform'
ONE_INC = 'one_big'
TWO_INC = 'two_big'
ONE_DEC = 'one_small'
TWO_DEC = 'two_small'

BALANCED = 'balanced'
HALPERN_BRUNO = 'halpern_bruno'

g_mode_to_color = {
        UNIFORM : 'black',
        ONE_INC : 'red',
        TWO_INC : 'orange',
        ONE_DEC : 'green',
        TWO_DEC : 'blue',
        }

g_ordered_modes = (
        UNIFORM,
        ONE_INC,
        TWO_INC,
        ONE_DEC,
        TWO_DEC,
        )

TOPLEFT = 'topleft'
BOTTOMLEFT = 'bottomleft'
TOPRIGHT = 'topright'
BOTTOMRIGHT = 'bottomright'


def get_form():
    """
    @return: the body of a form
    """
    # define the form objects
    form_objects = [
            Form.Integer('nstates', 'number of states', 4, low=2, high=9),
            Form.Float('t_low', 'initial time',
                0, low_inclusive=0),
            Form.Float('t_high', 'final time',
                '0.8', low_exclusive=0),
            Form.Float('sel_surr',
                'unitless deviation from uniformity (between 0 and 1)',
                '0.95', low_inclusive=0, high_exclusive=1),
            Form.RadioGroup('selection', 'selection approximation', [
                Form.RadioItem(BALANCED, 'WAG-like f=1/2'),
                Form.RadioItem(HALPERN_BRUNO, 'Halpern-Bruno', True)]),
            Form.CheckGroup('distribution', 'stationary distribution', [
                Form.CheckItem(UNIFORM, UNIFORM, True),
                Form.CheckItem(ONE_INC, ONE_INC.replace('_', ' '), True),
                Form.CheckItem(TWO_INC, TWO_INC.replace('_', ' '), True),
                Form.CheckItem(ONE_DEC, ONE_DEC.replace('_', ' '), True),
                Form.CheckItem(TWO_DEC, TWO_DEC.replace('_', ' '))]),
            Form.RadioGroup('legend_placement', 'plot legend location', [
                Form.RadioItem(TOPLEFT, 'top left'),
                Form.RadioItem(BOTTOMLEFT, 'bottom left'),
                Form.RadioItem(TOPRIGHT, 'top right', True),
                Form.RadioItem(BOTTOMRIGHT, 'bottom right')]),
            Form.ImageFormat(),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.Image('plot')

def get_distn_uniform(n, t):
    v = np.ones(n, dtype=float)
    return v / np.sum(v)

def get_distn_one_inc(n, t):
    v = (1-t) * np.ones(n, dtype=float)
    v[0] = 1
    return v / np.sum(v)

def get_distn_two_inc(n, t):
    v = (1-t) * np.ones(n, dtype=float)
    v[0] = 1
    v[1] = 1
    return v / np.sum(v)

def get_distn_one_dec(n, t):
    v = np.ones(n, dtype=float)
    v[0] -= t
    return v / np.sum(v)

def get_distn_two_dec(n, t):
    v = np.ones(n, dtype=float)
    v[0] -= t
    v[1] -= t
    return v / np.sum(v)

def make_table(args, distn_modes):
    """
    Make outputs to pass to RUtil.get_table_string.
    @param args: user args
    @param distn_modes: ordered distribution modes
    @return: matrix, headers
    """
    # define some variables
    t_low = args.t_low
    t_high = args.t_high
    if t_high <= t_low:
        raise ValueError('low time must be smaller than high time')
    ntimes = 100
    incr = (t_high - t_low) / (ntimes - 1)
    n = args.nstates
    # define some tables
    distn_mode_to_f = {
            UNIFORM : get_distn_uniform,
            ONE_INC : get_distn_one_inc,
            TWO_INC : get_distn_two_inc,
            ONE_DEC : get_distn_one_dec,
            TWO_DEC : get_distn_two_dec}
    selection_mode_to_f = {
            BALANCED : to_gtr_balanced,
            HALPERN_BRUNO : to_gtr_halpern_bruno}
    # define the selection modes and calculators
    selection_f = selection_mode_to_f[args.selection]
    distn_fs = [distn_mode_to_f[m] for m in distn_modes]
    # define the headers
    headers = ['t'] + [s.replace('_', '.') for s in distn_modes]
    # define the numbers in the table
    S = np.ones((n, n), dtype=float)
    S -= np.diag(np.sum(S, axis=1))
    arr = []
    for i in range(ntimes):
        t = t_low + i * incr
        row = [t]
        for distn_f in distn_fs:
            v = distn_f(n, args.sel_surr)
            R = selection_f(S, v)
            expected_log_ll_ratio = get_expected_ll_ratio(R, t)
            row.append(expected_log_ll_ratio)
        arr.append(row)
    return np.array(arr), headers

def to_gtr_balanced(S, v):
    """
    @param S: symmetric rate matrix
    @param v: target stationary distribution
    @return: a rate matrix with the target stationary distribution
    """
    # get the number of states
    n = len(v)
    # copy the symmetric rate matrix
    R = S.copy()
    # adjust the entries of the rate matrix
    for i in range(n):
        for j in range(n):
            R[i, j] *= math.sqrt(v[j] / v[i])
    # reset the diagonal entries of the rate matrix
    R -= np.diag(np.sum(R, axis=1))
    return R

def to_gtr_halpern_bruno(S, v):
    """
    @param S: symmetric rate matrix
    @param v: target stationary distribution
    @return: a rate matrix with the target stationary distribution
    """
    p = mrate.R_to_distn(S)
    # get the number of states
    n = len(v)
    # copy the symmetric rate matrix
    R = S.copy()
    # adjust the entries of the rate matrix
    for a in range(n):
        for b in range(n):
            if a != b:
                # This equation is unnecessarily verbose
                # due to symmetry of S.
                # It should also work for asymmetric input rate matrices.
                tau = (v[b] / p[b]) / (v[a] / p[a])
                if not np.allclose(tau, 1):
                    R[a, b] *= math.log(tau) / (1 - 1/tau)
    # reset the diagonal entries of the rate matrix
    R -= np.diag(np.sum(R, axis=1))
    return R

def get_expected_ll_ratio(R, t):
    # define the number of states
    n = len(R)
    # define the transition matrix
    P = scipy.linalg.expm(R*t)
    # define the stationary distribution
    p = mrate.R_to_distn(R)
    # get the expected log likelihood ratio
    accum = 0
    for i in range(n):
        for j in range(n):
            if p[i] and P[i, j]:
                coeff = p[i] * P[i, j]
                # cancel the p[i] in the numerator and denominator
                #numerator = p[i] * P[i, j]
                #denominator = p[i] * p[j]
                numerator = P[i, j]
                denominator = p[j]
                value = coeff * math.log(numerator / denominator)
                if not np.allclose(np.imag(value), 0):
                    raise ValueError('rogue imaginary number')
                accum += np.real(value)
    return accum

def get_response_content(fs):
    distn_modes = [x for x in g_ordered_modes if x in fs.distribution]
    if not distn_modes:
        raise ValueError('no distribution mode was specified')
    colors = [g_mode_to_color[m] for m in distn_modes]
    arr, headers = make_table(fs, distn_modes)
    distn_headers = headers[1:]
    # Get the largest value in the array,
    # skipping the first column.
    arrmax = np.max(arr[:,1:])
    # write the R script body
    out = StringIO()
    ylim = RUtil.mk_call_str('c', 0, arrmax + 0.1)
    sel_str = {
            BALANCED : 'f=1/2',
            HALPERN_BRUNO : 'Halpern-Bruno',
            }[fs.selection]
    print >> out, RUtil.mk_call_str(
            'plot',
            'my.table$t',
            'my.table$%s' % distn_headers[0],
            type='"n"',
            ylim=ylim,
            xlab='"time"',
            ylab='"expected log L-ratio"',
            main='"Effect of selection (%s) on log L-ratio for %d states"' % (sel_str, fs.nstates),
            )
    for c, header in zip(colors, distn_headers):
        print >> out, RUtil.mk_call_str(
                'lines',
                'my.table$t',
                'my.table$%s' % header,
                col='"%s"' % c,
                )
    mode_names = [s.replace('_', ' ') for s in distn_modes]
    legend_name_str = 'c(' + ', '.join('"%s"' % s for s in mode_names) + ')'
    legend_col_str = 'c(' + ', '.join('"%s"' % s for s in colors) + ')'
    legend_lty_str = 'c(' + ', '.join(['1']*len(distn_modes)) + ')'
    print >> out, RUtil.mk_call_str(
            'legend',
            '"%s"' % fs.legend_placement,
            legend_name_str,
            col=legend_col_str,
            lty=legend_lty_str,
            )
    script_body = out.getvalue()
    # create the R plot image
    table_string = RUtil.get_table_string(arr, headers)
    device_name = Form.g_imageformat_to_r_function[fs.imageformat]
    retcode, r_out, r_err, image_data = RUtil.run_plotter(
            table_string, script_body, device_name)
    if retcode:
        raise RUtil.RError(r_err)
    return image_data

