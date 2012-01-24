"""
Plot statistics of mutation-selection processes with random selection.

Look at histograms of ratios of statistics of
random mutation-selection balance rate matrices to
those of the pure mutation rate matrix.
In particular consider the expected rate (ER),
the normalized saturation rate (NSR),
and their product.
"""

from StringIO import StringIO
import argparse
import math
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

# This script was originally written
# as a bunch of separate R tikz plots that were pasted together
# into a latex file.
# But it has been rewritten so that all of the plots
# are generated as a single R image which
# can be converted in a better way to other image formats.

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
                100, low=1, high=10000),
            Form.RadioGroup('sel_var', 'selection parameter variance', [
                Form.RadioItem('low_var', 'low variance'),
                Form.RadioItem('medium_var', 'medium variance', True),
                Form.RadioItem('high_var', 'high variance'),
                Form.RadioItem('really_high_var', 'really high variance')]),
            Form.RadioGroup('sel_skew', 'selection parameter skew', [
                Form.RadioItem('neg_skew', 'some are very fit'),
                Form.RadioItem('no_skew', 'no skew', True),
                Form.RadioItem('pos_skew', 'some are very unfit')]),
            #Form.RadioGroup('legend_placement', 'plot legend location', [
                #Form.RadioItem(TOPLEFT, 'top left'),
                #Form.RadioItem(BOTTOMLEFT, 'bottom left'),
                #Form.RadioItem(TOPRIGHT, 'top right', True),
                #Form.RadioItem(BOTTOMRIGHT, 'bottom right')]),
            #Form.LatexFormat(),
            Form.ImageFormat(),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    #return FormOut.Latex('random-selection-report')
    return FormOut.Image('random-selection-plot')

def get_r_tikz_stub():
    user_script = RUtil.g_stub
    device_name = 'tikz'
    retcode, r_out, r_err, tikz_code = RUtil.run_plotter_no_table(
            user_script, device_name)
    if retcode:
        raise RUtil.RError(r_err)
    return tikz_code

def get_statistic_ratios(Q_mut, Q_sels):
    """
    @param Q_mut: mutation rate matrix
    @param Q_sels: mutations-selection balance rate matrices
    @return: ER_ratios, NSR_ratios, ER_NSR_ratios
    """
    ER_mut = mrate.R_to_total_rate(Q_mut)
    ER_sels = [mrate.R_to_total_rate(Q) for Q in Q_sels]
    ER_ratios = [ER_sel / ER_mut for ER_sel in ER_sels]
    ER_NSR_mut = 1 / mrate.R_to_relaxation_time(Q_mut)
    ER_NSR_sels = [1 / mrate.R_to_relaxation_time(Q) for Q in Q_sels]
    ER_NSR_ratios = [ER_NSR_sel / ER_NSR_mut for ER_NSR_sel in ER_NSR_sels]
    NSR_ratios = [a / b for a, b in zip(ER_NSR_ratios, ER_ratios)]
    # do some extra investigation
    """
    nsels = len(Q_sels)
    for i in range(nsels):
        if ER_NSR_ratios[i] < 1:
            print 'found a slower-decaying mutation-selection matrix:'
            print Q_sels[i]
            print
    print
    print 'ER_mut:'
    print ER_mut
    print
    print 'ER_NSR_mut:'
    print ER_NSR_mut
    print
    print 'ER_sels:'
    for x in ER_sels:
        print x
    print
    print 'ER_NSR_sels:'
    for x in ER_NSR_sels:
        print x
    print
    """
    return ER_ratios, NSR_ratios, ER_NSR_ratios

def get_r_tikz_script(nsels, name):
    """
    This is obsolete because I am now using pure R output.
    @param nsels: the number of mutation-selection balance matrices
    @param name: the name of the variable whose histogram should be plotted
    @return: tikz code corresponding to an R plot
    """
    out = StringIO()
    print >> out, RUtil.mk_call_str(
            'hist',
            'my.table$%s' % name,
            xlab = '"%s"' % name,
            ylab = '"counts"',
            main='"%s; N=%d"' % (name, nsels))
    return out.getvalue()

def get_r_comboscript(nsels, names):
    """
    Do all of the plots at once.
    @param nsels: the number of mutation-selection balance matrices
    @param names: the names of the variable whose histogram should be plotted
    @return: tikz code corresponding to an R plot
    """
    out = StringIO()
    print >> out, 'par(mfrow=c(3,1))'
    for name in names:
        print >> out, RUtil.mk_call_str(
                'hist',
                'my.table$%s' % name,
                xlab = '"%s"' % name,
                ylab = '"counts"',
                main='"%s; N=%d"' % (name, nsels))
    return out.getvalue()

def get_qmut_qsels(fs):
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
            v = 1.0
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
    return Q_mut, Q_sels

def get_latex_documentbody(fs):
    """
    This is obsolete because I am now using pure R output.
    The latex documentbody should have a bunch of tikz pieces in it.
    Each tikz piece should have been generated from R.
    """
    Q_mut, Q_sels = get_qmut_qsels(fs)
    # compute the statistics
    ER_ratios, NSR_ratios, ER_NSR_ratios  = get_statistic_ratios(Q_mut, Q_sels)
    M = zip(*(ER_ratios, NSR_ratios, ER_NSR_ratios))
    column_headers = ('ER.ratio', 'NSR.ratio', 'ER.times.NSR.ratio')
    table_string = RUtil.get_table_string(M, column_headers)
    nsels = len(Q_sels)
    # define the R scripts
    scripts = []
    for name in column_headers:
        scripts.append(get_r_tikz_script(nsels, name))
    # get the tikz codes from R, for each histogram
    retcode, r_out, r_err, tikz_code_list = RUtil.run_plotter_multiple_scripts(
            table_string, scripts, 'tikz',
            width=3, height=2)
    if retcode:
        raise RUtil.RError(r_err)
    #
    # show some timings
    print 'R did not fail, but here is its stderr:'
    print r_err
    #
    # write the latex code
    out = StringIO()
    #print >> out, '\\pagestyle{empty}'
    for tikz_code in tikz_code_list:
        print >> out, tikz_code
    # return the latex code, consisting mainly of a bunch of tikz plots
    return out.getvalue()

def get_response_content_latex(fs):
    """
    This is obsolete because I am now using pure R output.
    """
    requested_documentclass = 'standalone'
    document_body = get_latex_documentbody(fs)
    latexformat = fs.latexformat
    packages = ('tikz', 'verbatim')
    preamble = ''
    return latexutil.get_response(
            requested_documentclass, document_body, latexformat,
            packages, preamble)

def get_response_content(fs):
    Q_mut, Q_sels = get_qmut_qsels(fs)
    # compute the statistics
    ER_ratios, NSR_ratios, ER_NSR_ratios  = get_statistic_ratios(Q_mut, Q_sels)
    M = zip(*(ER_ratios, NSR_ratios, ER_NSR_ratios))
    column_headers = ('ER.ratio', 'NSR.ratio', 'ER.times.NSR.ratio')
    table_string = RUtil.get_table_string(M, column_headers)
    nsels = len(Q_sels)
    # get the R script
    comboscript = get_r_comboscript(nsels, column_headers)
    # create the R plot image 
    device_name = Form.g_imageformat_to_r_function[fs.imageformat] 
    retcode, r_out, r_err, image_data = RUtil.run_plotter( 
        table_string, comboscript, device_name) 
    if retcode: 
        raise RUtil.RError(r_err) 
    return image_data 

