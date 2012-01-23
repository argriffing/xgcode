"""
Plot mutual information over time and over many mutuation-selection processes.
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
import ctmcmi
import latexutil
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

def sample_low_variance_params(n):
    variance = 0.01
    mu = 0.0
    s = math.sqrt(variance)
    return [random.gauss(mu, s) for i in range(n)]

def sample_high_variance_params(n):
    variance = 0.1
    mu = 0.0
    s = math.sqrt(variance)
    return [random.gauss(mu, s) for i in range(n)]

def sample_neg_skewed_params(n):
    variance = 0.02
    mu = 0.0
    s = math.sqr(variance)
    rate = s
    return [random.gauss(mu, s) - random.expovariate(rate) for i in range(n)]

def sample_pos_skewed_params(n):
    variance = 0.02
    mu = 0.0
    s = math.sqr(variance)
    rate = s
    return [random.gauss(mu, s) + random.expovariate(rate) for i in range(n)]

def get_dense_sequence_rate_matrix(nsites):
    """
    Each sequences changes to each other sequence at the same rate.
    The matrix is normalized by expected rate.
    """
    dna = 4
    nstates = dna**nsites
    R = np.ones((nstates, nstates))
    for i in range(nstates):
        R[i, i] = -(nstates - 1)
    uniform_pi = np.reciprocal(nstates * np.ones(nstates))
    expected_rate = -sum(uniform_pi[i] * R[i, i] for i in range(nstates))
    return R / expected_rate

def get_sparse_sequence_rate_matrix(nsites):
    """
    Sites change change independently.
    The matrix is normalized by expected rate.
    """
    dna = 4
    nstates = dna**nsites
    R = np.zeros((nstates, nstates))
    for alpha in itertools.product(range(dna), repeat=nsites):
        for beta in itertools.product(range(dna), repeat=nsites):
            alpha_index = [alpha[i]*(dna ** i) for i in range(nsites)]
            beta_index = [beta[i]*(dna ** i) for i in range(nsites)]
            hamming_dist = sum(1 for a, b in zip(alpha, beta) if a != b)
            if hamming_dist == 1:
                R[alpha_index, beta_index] = -1
    for i in range(n):
        R[i, i] = -np.sum(R[i])
    uniform_pi = np.reciprocal(nstates * np.ones(nstates))
    expected_rate = -sum(uniform_pi[i] * R[i, i] for i in range(nstates))
    return R / expected_rate



def get_form():
    """
    @return: the body of a form
    """
    # define the form objects
    form_objects = [
            Form.Integer('nsites', 'number of nucleotide sites',
                3, low=1, high=3),
            Form.RadioGroup('mut_type', 'mutation matrix form', [
                Form.RadioItem('sparse', 'J-C one site at a time', True),
                Form.RadioItem('dense', 'all seq-seq rates equal')]),
            Form.Integer('nselections', 'number of mutation-selection samples',
                100, low=3, high=1000),
            Form.RadioGroup('sel_var', 'selection parameter variance', [
                Form.RadioItem('low_var', 'low variance'),
                Form.RadioItem('medium_var', 'medium variance', True),
                Form.RadioItem('high_var', 'high variance')]),
            Form.RadioGroup('sel_skew', 'selection parameter skew', [
                Form.RadioItem('neg_skew', 'some very fit'),
                Form.RadioItem('no_skew', 'no skew', True),
                Form.RadioItem('pos_skew', 'some very unfit')]),
            Form.Float('t_low', 'initial time',
                '0.001', low_exclusive=0),
            Form.Float('t_high', 'final time',
                '3.0', low_exclusive=0),
            Form.Integer('ntimes', 'sample this many time points',
                100, low=3, high=100),
            Form.CheckGroup('surrogate_functions',
                'show mutual information correlation with these functions', [
                    Form.CheckItem('mi_diag_approx',
                        'diagonal contribution to mutual info (approx)', True),
                    Form.CheckItem('mi_diag',
                        'diagonal contribution to mutual info', True),
                    Form.CheckItem('large_t_approx',
                        'large time approximation to mutual info', True),
                    Form.CheckItem('eigenvalue',
                        'exponential function of second eigenvalue', True),
                    Form.CheckItem('expected_rate',
                        'exponential function of expected rate', True)]),
            #Form.RadioGroup('legend_placement', 'plot legend location', [
                #Form.RadioItem(TOPLEFT, 'top left'),
                #Form.RadioItem(BOTTOMLEFT, 'bottom left'),
                #Form.RadioItem(TOPRIGHT, 'top right', True),
                #Form.RadioItem(BOTTOMRIGHT, 'bottom right')]),
            Form.LatexFormat(),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.Latex('mutual-information-report')

def get_response_content(fs):
    requested_documentclass = 'article'
    document_body = 'o hai'
    latexformat = fs.latexformat
    packages = ()
    preamble = ''
    return latexutil.get_response(
            requested_documentclass, document_body, latexformat,
            packages, preamble)

def OBSOLETE_make_table(args, distn_modes):
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
        raise ValueError('invalid time point range')
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

def OBSOLETE_get_response_content(fs):
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

