"""
Plot correlations of functions related to mutual information.

The red curve is an approximation of the
diagonal contribution to the mutual information.
The orange curve is the exact
diagonal contribution to the mutual information.
The green curve is a good large-time approximation.
The blue curve is an exponential function of the lambda_2 eigenvalue.
The black curve is an exponential function of the expected rate.
The idea is that some of these functions could have a more intuitive meaning
or a closer connection to the spectral or parametric representations
of the rate matrices while being asymptotically related
to the mutual information.
"""

from StringIO import StringIO
import argparse
import math
import itertools
import random

import numpy as np
import scipy
from scipy import linalg
from scipy import stats

import Form
import FormOut
from MatrixUtil import ndot
import mrate
import ctmcmi
import latexutil
import tikz
import RUtil

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
            Form.Integer('nresidues', 'number of states per site',
                2, low=2, high=256),
            Form.Integer('nsites', 'number of sites',
                2, low=1, high=8),
            Form.Integer('nselections', 'number of mutation-selection samples',
                100, low=100, high=1000),
            Form.RadioGroup('sel_var', 'selection parameter variance', [
                Form.RadioItem('low_var', 'low variance'),
                Form.RadioItem('medium_var', 'medium variance', True),
                Form.RadioItem('high_var', 'high variance'),
                Form.RadioItem('really_high_var', 'really high variance')]),
            Form.RadioGroup('sel_skew', 'selection parameter skew', [
                Form.RadioItem('neg_skew', 'some are very fit'),
                Form.RadioItem('no_skew', 'no skew', True),
                Form.RadioItem('pos_skew', 'some are very unfit')]),
            Form.Integer('ntimes', 'sample this many time points',
                10, low=3, high=100),
            Form.FloatInterval(
                't_low', 't_high', 'divtime interval',
                '0.001', '4.0', low_exclusive=0),
            #Form.CheckGroup('surrogate_functions',
                #'show mutual information correlation with these functions', [
                    #Form.CheckItem('mi_diag_approx',
                        #'diagonal contribution to mutual info (approx)', True),
                    #Form.CheckItem('mi_diag',
                        #'diagonal contribution to mutual info', True),
                    #Form.CheckItem('large_t_approx',
                        #'large time approximation to mutual info', True),
                    #Form.CheckItem('eigenvalue',
                        #'exponential function of second eigenvalue', True),
                    #Form.CheckItem('expected_rate',
                        #'exponential function of expected rate', True)]),
            #Form.RadioGroup('legend_placement', 'plot legend location', [
                #Form.RadioItem(TOPLEFT, 'top left'),
                #Form.RadioItem(BOTTOMLEFT, 'bottom left'),
                #Form.RadioItem(TOPRIGHT, 'top right', True),
                #Form.RadioItem(BOTTOMRIGHT, 'bottom right')]),
            #Form.LatexFormat(),
            Form.ImageFormat()]
    return form_objects

def get_form_out():
    #return FormOut.Latex('mutual-information-report')
    return FormOut.Image('mutual-information-correlates')

def get_r_tikz_stub():
    user_script = RUtil.g_stub
    device_name = 'tikz'
    retcode, r_out, r_err, tikz_code = RUtil.run_plotter_no_table(
            user_script, device_name)
    if retcode:
        raise RUtil.RError(r_err)
    return tikz_code

def get_time_point_summary(Q_mut, Q_sels, t):
    """
    @param Q_mut: the mutation rate matrix
    @param Q_sels: sequence of mutation-selection rate matrices
    @param t: the time point under consideration
    @return: a sequence of statistics
    """
    # Compute the following statistics at this time point:
    # t
    # mutation MI
    # selection MI max
    # selection MI high
    # selection MI mean
    # selection MI low
    # selection MI min
    # correlation fn 1
    # correlation fn 2
    # correlation fn 3
    # correlation fn 4
    # correlation fn 5
    # proportion sign agreement fn 1
    # proportion sign agreement fn 2
    # proportion sign agreement fn 3
    # proportion sign agreement fn 4
    # proportion sign agreement fn 5
    # informativeness fn 1
    # informativeness fn 2
    # informativeness fn 3
    # informativeness fn 4
    # informativeness fn 5
    #
    # First compute the mutual information for mut and mut-sel.
    nsels = len(Q_sels)
    mi_mut = ctmcmi.get_mutual_information(Q_mut, t)
    mi_sels = [ctmcmi.get_mutual_information(Q, t) for Q in Q_sels]
    mi_signs = [1 if mi_sel > mi_mut else -1 for mi_sel in mi_sels]
    # Now compute some other functions
    v0 = [ctmcmi.get_mutual_information_small_approx_c(Q, t) for Q in Q_sels]
    v1 = [ctmcmi.get_mutual_information_small_approx(Q, t) for Q in Q_sels]
    v2 = [ctmcmi.get_mutual_information_approx_c(Q, t) for Q in Q_sels]
    v3 = [math.exp(-2*t/mrate.R_to_relaxation_time(Q)) for Q in Q_sels]
    v4 = [math.exp(-t*mrate.Q_to_expected_rate(Q)) for Q in Q_sels]
    # Now that we have computed all of the vectors at this time point,
    # we can compute the statistics that we want to report.
    statistics = []
    statistics.append(t)
    statistics.append(mi_mut)
    # add the mutual information statistics
    sorted_mi = sorted(mi_sels)
    n_extreme = nsels / 20
    statistics.append(sorted_mi[-1])
    statistics.append(sorted_mi[-n_extreme])
    statistics.append(sum(sorted_mi) / nsels)
    statistics.append(sorted_mi[n_extreme-1])
    statistics.append(sorted_mi[0])
    # add the correlations
    for v in (v0, v1, v2, v3, v4):
        r, p = scipy.stats.stats.pearsonr(v, mi_sels)
        statistics.append(r)
    # add the sign proportions
    for v in (v0, v1, v2, v3, v4):
        v_signs = [1 if value > mi_mut else -1 for value in v]
        total = sum(1 for a, b in zip(mi_signs, v_signs) if a == b)
        p = float(total) / nsels
        statistics.append(p)
    # add the informativenesses
    for v in (v0, v1, v2, v3, v4):
        v_signs = [1 if value > mi_mut else -1 for value in v]
        informativeness = 0
        for pair in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
            v_value, m_value = pair
            v_marginal_count = sum(1 for x in v_signs if x == v_value)
            m_marginal_count = sum(1 for x in mi_signs if x == m_value)
            joint_count = sum(1 for x in zip(v_signs, mi_signs) if x == pair)
            if joint_count:
                joint_prob = joint_count / float(nsels)
                a = math.log(joint_prob)
                b = math.log(v_marginal_count / float(nsels))
                c = math.log(m_marginal_count / float(nsels))
                informativeness += joint_prob * (a - b - c)
        statistics.append(informativeness)
    # return the statistics
    return statistics

g_time_stats_headers = (
        't', 'mut', 'mut.sel.max', 'mut.sel.high',
        'mut.sel.mean', 'mut.sel.low', 'mut.sel.min',
        'corr.mi.diag.approx', 'corr.mi.diag', 'corr.large.t.approx',
        'corr.expon.eigen', 'corr.expon.e.rate',
        'prop.mi.diag.approx', 'prop.mi.diag', 'prop.large.t.approx',
        'prop.expon.eigen', 'prop.expon.e.rate',
        'info.mi.diag.approx', 'info.mi.diag', 'info.large.t.approx',
        'info.expon.eigen', 'info.expon.e.rate')

def get_r_tikz_mi_plot_script(nsels, time_stats):
    """
    At each time point plot mutual information for all matrices.
    @param time_stats: a list of stats for each time point
    @return: tikz code corresponding to an R plot
    """
    out = StringIO()
    time_stats_trans = zip(*time_stats)
    mi_mut = time_stats_trans[1]
    mi_min_sels = time_stats_trans[6]
    mi_max_sels = time_stats_trans[2]
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
            main='"MI for mut process and %d mut.sel processes"' % nsels)
    colors = ('red', 'blue', 'green', 'black', 'green', 'blue')
    plot_indices = (1, 2, 3, 4, 5, 6)
    for c, plot_index in zip(colors, plot_indices):
        header = g_time_stats_headers[plot_index]
        print >> out, RUtil.mk_call_str(
                'lines',
                'my.table$t',
                'my.table$%s' % header,
                col='"%s"' % c)
    return out.getvalue()

def get_r_tikz_corr_plot(nsels, time_stats):
    """
    @param time_stats: a list of stats for each time point
    @return: tikz code corresponding to an R plot
    """
    out = StringIO()
    time_stats_trans = zip(*time_stats)
    y_low = -1
    y_high = 1
    ylim = RUtil.mk_call_str('c', y_low, y_high)
    print >> out, RUtil.mk_call_str(
            'plot',
            'my.table$t',
            'my.table$corr.mi.diag.approx',
            type='"n"',
            ylim=ylim,
            xlab='"time"',
            ylab='"correlation"',
            main='"correlation with mutual information"')
    colors = ('red', 'orange', 'green', 'blue', 'black')
    plot_indices = (7, 8, 9, 10, 11)
    for c, plot_index in zip(colors, plot_indices):
        header = g_time_stats_headers[plot_index]
        print >> out, RUtil.mk_call_str(
                'lines',
                'my.table$t',
                'my.table$%s' % header,
                col='"%s"' % c)
    return out.getvalue()

def get_r_tikz_prop_plot(nsels, time_stats):
    """
    @param time_stats: a list of stats for each time point
    @return: tikz code corresponding to an R plot
    """
    out = StringIO()
    time_stats_trans = zip(*time_stats)
    y_low = 0
    y_high = 1
    ylim = RUtil.mk_call_str('c', y_low, y_high)
    print >> out, RUtil.mk_call_str(
            'plot',
            'my.table$t',
            'my.table$prop.mi.diag.approx',
            type='"n"',
            ylim=ylim,
            xlab='"time"',
            ylab='"proportion"',
            main='"proportion of same sign difference as MI"')
    colors = ('red', 'orange', 'green', 'blue', 'black')
    plot_indices = (12, 13, 14, 15, 16)
    for c, plot_index in zip(colors, plot_indices):
        header = g_time_stats_headers[plot_index]
        print >> out, RUtil.mk_call_str(
                'lines',
                'my.table$t',
                'my.table$%s' % header,
                col='"%s"' % c)
    return out.getvalue()

def get_r_tikz_info_plot(nsels, time_stats):
    """
    @param time_stats: a list of stats for each time point
    @return: tikz code corresponding to an R plot
    """
    out = StringIO()
    time_stats_trans = zip(*time_stats)
    y_low = 0
    y_high = math.log(2)
    ylim = RUtil.mk_call_str('c', y_low, y_high)
    print >> out, RUtil.mk_call_str(
            'plot',
            'my.table$t',
            'my.table$info.mi.diag.approx',
            type='"n"',
            ylim=ylim,
            xlab='"time"',
            ylab='"info"',
            main='"informativeness with respect to MI"')
    colors = ('red', 'orange', 'green', 'blue', 'black')
    plot_indices = (17, 18, 19, 20, 21)
    for c, plot_index in zip(colors, plot_indices):
        header = g_time_stats_headers[plot_index]
        print >> out, RUtil.mk_call_str(
                'lines',
                'my.table$t',
                'my.table$%s' % header,
                col='"%s"' % c)
    return out.getvalue()

def get_table_string_and_scripts(fs):
    """
    The latex documentbody should have a bunch of tikz pieces in it.
    Each tikz piece should have been generated from R.
    """
    nstates = fs.nresidues ** fs.nsites
    if nstates > 256:
        raise ValueError('the mutation rate matrix is too big')
    # get the mutation matrix
    Q_mut = mrate.get_sparse_sequence_rate_matrix(fs.nresidues, fs.nsites)
    # sample a bunch of mutation-selection rate matrices
    Q_sels = []
    for selection_index in range(fs.nselections):
        # sample the selection parameters
        if fs.low_var:
            v = 0.2
        elif fs.medium_var:
            v = 1
        elif fs.high_var:
            v = 5.0
        elif fs.really_high_var:
            v = 25.0
        s = math.sqrt(v)
        if fs.neg_skew:
            sels = [-random.expovariate(1/s) for i in range(nstates)]
        elif fs.no_skew:
            sels = [random.gauss(0, s) for i in range(nstates)]
        elif fs.pos_skew:
            sels = [random.expovariate(1/s) for i in range(nstates)]
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
    # compute the statistics
    nsels = len(Q_sels)
    time_stats = [get_time_point_summary(Q_mut, Q_sels, t) for t in times]
    # get the R scripts
    scripts = [
            #get_r_tikz_mi_plot(nsels, time_stats),
            get_r_tikz_corr_plot(nsels, time_stats),
            get_r_tikz_prop_plot(nsels, time_stats),
            get_r_tikz_info_plot(nsels, time_stats)]
    table_string = RUtil.get_table_string(time_stats, g_time_stats_headers)
    return table_string, scripts

def get_latex_documentbody(fs):
    """
    This is obsolete.
    """
    out = StringIO()
    table_string, scripts = get_table_string_and_scripts(fs)
    for script in scripts:
        retcode, r_out, r_err, tikz_code = RUtil.run_plotter(
                table_string, script, 'tikz',
                width=5, height=5)
        if retcode:
            raise RUtil.RError(r_err)
        print >> out, tikz_code
    return out.getvalue()

def get_response_content_latex(fs):
    requested_documentclass = 'article'
    document_body = get_latex_documentbody(fs)
    latexformat = fs.latexformat
    packages = ('tikz', 'verbatim')
    preamble = ''
    return latexutil.get_response(
            requested_documentclass, document_body, latexformat,
            packages, preamble)

def get_response_content(fs):
    # get the table string and scripts
    table_string, scripts = get_table_string_and_scripts(fs)
    # create a comboscript
    out = StringIO()
    print >> out, 'par(mfrow=c(3,1))'
    for script in scripts:
        print >> out, script
    comboscript = out.getvalue()
    # create the R plot image 
    device_name = Form.g_imageformat_to_r_function[fs.imageformat] 
    retcode, r_out, r_err, image_data = RUtil.run_plotter( 
        table_string, comboscript, device_name) 
    if retcode: 
        raise RUtil.RError(r_err) 
    return image_data 

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
            BALANCED : mrate.to_gtr_balanced,
            HALPERN_BRUNO : mrate.to_gtr_halpern_bruno}
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

