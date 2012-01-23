"""
Plot mutual information over time and over many mutuation-selection processes.
"""

from StringIO import StringIO
import argparse
import math
import itertools
import random

import numpy as np
import scipy
from scipy import linalg

import Form
import FormOut
from MatrixUtil import ndot
import mrate
import ctmcmi
import latexutil
import tikz
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
            alpha_index = sum(alpha[i]*(dna ** i) for i in range(nsites))
            beta_index = sum(beta[i]*(dna ** i) for i in range(nsites))
            hamming_dist = sum(1 for a, b in zip(alpha, beta) if a != b)
            if hamming_dist == 1:
                R[alpha_index, beta_index] = 1
    for i in range(nstates):
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
                100, low=100, high=1000),
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

def get_r_tikz_stub():
    user_script = RUtil.g_stub
    device_name = 'tikz'
    retcode, r_out, r_err, tikz_code = RUtil.run_plotter_no_table(
            user_script, device_name)
    if retcode:
        raise RUtil.RError(r_err)
    return tikz_code

def get_r_tikz_mi_plot(times, Q_mut, Q_sels):
    """
    At each time point plot mutual information for all matrices.
    @param times: a list of times
    @param Q_mut: a mutation rate matrix
    @param Q_sels: a sequence of mutation-selection rate matrices
    @return: tikz code corresponding to an R plot
    """
    mi_mut = [ctmcmi.get_mutual_information(Q_mut, t) for t in times]
    mi_max_sels = []
    mi_high_sels = []
    mi_mean_sels = []
    mi_low_sels = []
    mi_min_sels = []
    n_extreme = len(Q_sels) / 20
    for t in times:
        data = sorted(ctmcmi.get_mutual_information(Q, t) for Q in Q_sels)
        mi_max_sels.append(data[-1])
        mi_high_sels.append(data[-n_extreme])
        mi_mean_sels.append(sum(data) / len(Q_sels))
        mi_low_sels.append(data[n_extreme-1])
        mi_min_sels.append(data[0])
    # define the R data table
    column_headers = ('t', 'mut', 'mut.sel.max', 'mut.sel.high',
            'mut.sel.mean', 'mut.sel.low', 'mut.sel.min')
    M = zip(*(times, mi_mut, mi_max_sels, mi_high_sels,
        mi_mean_sels, mi_low_sels, mi_min_sels))
    # define the script content
    out = StringIO()
    y_low = min(mi_min_sels + mi_mut)
    y_high = max(mi_max_sels + mi_mut)
    ylim = RUtil.mk_call_str('c', y_low, y_high)
    print >> out, RUtil.mk_call_str(
            'plot',
            'my.table$t',
            'my.table$mut',
            type='"n"',
            ylim=ylim,
            xlab='"time"',
            ylab='"MI"',
            main='"MI for mut process and %d mut.sel processes"' % len(Q_sels))
    colors = ('red', 'blue', 'green', 'black', 'green', 'blue')
    for c, header in zip(colors, column_headers[1:]):
        print >> out, RUtil.mk_call_str(
                'lines',
                'my.table$t',
                'my.table$%s' % header,
                col='"%s"' % c)
    user_script_content = out.getvalue()
    # call R to get the tikz code
    table_string = RUtil.get_table_string(M, column_headers)
    retcode, r_out, r_err, tikz_code = RUtil.run_plotter(
            table_string, user_script_content, 'tikz')
    if retcode:
        raise RUtil.RError(r_err)
    return tikz_code

def get_latex_documentbody(fs):
    """
    The latex documentbody should have a bunch of tikz pieces in it.
    Each tikz piece should have been generated from R.
    """
    nstates = 4 ** fs.nsites
    # get the mutation matrix
    if fs.sparse:
        Q_mut = get_sparse_sequence_rate_matrix(fs.nsites)
    elif fs.dense:
        Q_mut = get_dense_sequence_rate_matrix(fs.nsites)
    # sample a bunch of mutation-selection rate matrices
    Q_sels = []
    for selection_index in range(fs.nselections):
        # sample the selection parameters
        if fs.low_var:
            v = 0.1
        elif fs.medium_var:
            v = 0.4
        elif fs.high_var:
            v = 1.0
        s = math.sqrt(v)
        if fs.neg_skew:
            sels = [-random.expovariate(s) for i in range(nstates)]
        elif fs.no_skew:
            sels = [random.gauss(0, s) for i in range(nstates)]
        elif fs.pos_skew:
            sels = [random.expovariate(s) for i in range(nstates)]
        # define the mutation-selection rate matrix using Halpern-Bruno
        Q = np.zeros_like(Q_mut)
        for i in range(nstates):
            for j in range(nstates):
                if i != j:
                    tau = math.exp(-(sels[j] - sels[i]))
                    coeff = math.log(tau) / (1 - 1/tau)
                    Q[i, j] = Q_mut[i, j] * coeff
        for i in range(nstates):
            Q[i, i] = -np.sum(Q[i])
        Q_sels.append(Q)
    # define the time points
    incr = (fs.t_high - fs.t_low) / (fs.ntimes - 1)
    times = [fs.t_low + i*incr for i in range(fs.ntimes)]
    # write the latex code
    out = StringIO()
    print >> out, 'Hello world!'
    print >> out
    print >> out, get_r_tikz_mi_plot(times, Q_mut, Q_sels)
    return out.getvalue()

def get_response_content(fs):
    requested_documentclass = 'article'
    document_body = get_latex_documentbody(fs)
    latexformat = fs.latexformat
    packages = ('tikz', 'verbatim')
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

