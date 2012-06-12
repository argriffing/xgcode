"""
Plot summaries of the behavior of mutual information under random selection.
"""

from StringIO import StringIO
from collections import defaultdict
import math
import random

import numpy as np

import Form
import FormOut
import mrate
import ctmcmi
import RUtil
import iterutils

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
            #Form.RadioGroup('legend_placement', 'plot legend location', [
                #Form.RadioItem(TOPLEFT, 'top left'),
                #Form.RadioItem(BOTTOMLEFT, 'bottom left'),
                #Form.RadioItem(TOPRIGHT, 'top right', True),
                #Form.RadioItem(BOTTOMRIGHT, 'bottom right')]),
            Form.ImageFormat()]
    return form_objects

def get_form_out():
    return FormOut.Image('mutual-information-correlates')

def get_time_point_summary(Q_mut, Q_sels, t):
    """
    @param Q_mut: the mutation rate matrix
    @param Q_sels: sequence of mutation-selection rate matrices
    @param t: the time point under consideration
    @return: a list of signs, and a sequence of statistics
    """
    # Compute the following statistics at this time point:
    # t
    # mutation MI
    # selection MI max
    # selection MI high
    # selection MI mean
    # selection MI low
    # selection MI min
    # proportion
    #
    # First compute the mutual information for mut and mut-sel.
    nsels = len(Q_sels)
    mi_mut = ctmcmi.get_mutual_information(Q_mut, t)
    mi_sels = [ctmcmi.get_mutual_information(Q, t) for Q in Q_sels]
    mi_signs = [1 if mi_sel > mi_mut else -1 for mi_sel in mi_sels]
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
    # add the proportion
    statistics.append(sum(1 for x in mi_signs if x == 1) / float(nsels))
    # return the statistics
    return mi_signs, statistics

g_time_stats_headers = (
        't', 'mut', 'mut.sel.max', 'mut.sel.high',
        'mut.sel.mean', 'mut.sel.low', 'mut.sel.min',
        'prop.sel.vs.mut')

def get_r_band_script(nsels, time_stats):
    """
    @param time_stats: a list of stats for each time point
    @return: R code
    """
    out = StringIO()
    time_stats_trans = zip(*time_stats)
    mi_mut = time_stats_trans[1]
    # set up the correctly sized plot
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
            main='"MI for mut process (red) and %d mut-sel processes"' % nsels)
    # draw a light gray polygon covering all selection mutual information
    print >> out, RUtil.mk_call_str(
            'polygon',
            'c(my.table$t, rev(my.table$t))',
            'c(my.table$mut.sel.max, rev(my.table$mut.sel.min))',
            col='"gray80"',
            border='NA')
    # draw a darker gray polygon covering most of selection mutual information
    print >> out, RUtil.mk_call_str(
            'polygon',
            'c(my.table$t, rev(my.table$t))',
            'c(my.table$mut.sel.high, rev(my.table$mut.sel.low))',
            col='"gray50"',
            border='NA')
    # draw the black line representing the mean selection mutual information
    print >> out, RUtil.mk_call_str(
            'lines',
            'my.table$t',
            'my.table$mut.sel.mean',
            col='"black"')
    # draw the red line representing the mutation mutual information
    print >> out, RUtil.mk_call_str(
            'lines',
            'my.table$t',
            'my.table$mut',
            col='"red"')
    return out.getvalue()

def get_r_prop_script(nsels, time_stats):
    """
    @param time_stats: a list of stats for each time point
    @return: R code
    """
    out = StringIO()
    time_stats_trans = zip(*time_stats)
    y_low = 0
    y_high = 1
    ylim = RUtil.mk_call_str('c', y_low, y_high)
    print >> out, RUtil.mk_call_str(
            'plot',
            'my.table$t',
            'my.table$prop.sel.vs.mut',
            type='"l"',
            ylim=ylim,
            xlab='"time"',
            ylab='"proportion"',
            main='"proportion of mut-sel MI greater than mutation MI"')
    return out.getvalue()

def get_r_cross_script(ncrossing_list):
    """
    @param time_stats: a list of stats for each time point
    @return: R code
    """
    out = StringIO()
    low = min(ncrossing_list)
    high = max(ncrossing_list)
    n_to_count = defaultdict(int)
    for n in ncrossing_list:
        n_to_count[n] += 1
    counts = [n_to_count[i] for i in range(low, high+1)]
    s = ', '.join('"%s"' % i for i in range(low, high+1))
    print >> out, RUtil.mk_call_str(
            'barplot',
            'c(' + ', '.join(str(x) for x in counts) + ')',
            'names.arg=c(' + s + ')',
            xlab='"number of crossings"',
            ylab='"frequency"',
            main='"number of times mut-sel MI crosses mut MI"')
    return out.getvalue()

def get_table_string_and_scripts(fs):
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
    pairs = [get_time_point_summary(Q_mut, Q_sels, t) for t in times]
    mi_sign_lists, time_stats = zip(*pairs)
    ncrossing_list = []
    # look at how the signs change over time for each selection sample
    for signs in zip(*mi_sign_lists):
        count = 0
        for sign_a, sign_b in iterutils.pairwise(signs):
            if sign_a != sign_b:
                count += 1
        ncrossing_list.append(count)
    # get the R scripts
    scripts = [
            get_r_band_script(nsels, time_stats),
            get_r_prop_script(nsels, time_stats),
            get_r_cross_script(ncrossing_list)]
    table_string = RUtil.get_table_string(time_stats, g_time_stats_headers)
    return table_string, scripts

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

