"""
Plot statistics of mutation-selection processes with random selection.

Look at histograms of ratios of statistics of
random mutation-selection balance rate matrices to
those of the pure mutation rate matrix.
In particular consider the expected rate (ER)
and the normalized saturation rate (NSR).
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


def get_form():
    """
    @return: the body of a form
    """
    # define the form objects
    form_objects = [
            Form.Integer('nresidues', 'number of states per site',
                4, low=2, high=64),
            Form.Integer('nsites', 'number of sites',
                2, low=1, high=3),
            Form.RadioGroup('mut_type', 'mutation matrix form', [
                Form.RadioItem('sparse', 'one site changes at a time', True),
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
            #Form.RadioGroup('legend_placement', 'plot legend location', [
                #Form.RadioItem(TOPLEFT, 'top left'),
                #Form.RadioItem(BOTTOMLEFT, 'bottom left'),
                #Form.RadioItem(TOPRIGHT, 'top right', True),
                #Form.RadioItem(BOTTOMRIGHT, 'bottom right')]),
            Form.LatexFormat(),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.Latex('random-selection-report')

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
    #print
    #print 'ER_mut:'
    #print ER_mut
    #print
    #print 'ER_NSR_mut:'
    #print ER_NSR_mut
    #print
    #print 'ER_sels:'
    #for x in ER_sels:
        #print x
    #print
    #print 'ER_NSR_sels:'
    #for x in ER_NSR_sels:
        #print x
    #print
    return ER_ratios, NSR_ratios, ER_NSR_ratios

def get_r_tikz_hist(nsels, table_string, name):
    """
    @param table_string: the R table string
    @param name: the name of the variable whose histogram should be plotted
    @return: tikz code corresponding to an R plot
    """
    # define the script content
    out = StringIO()
    print >> out, RUtil.mk_call_str(
            'hist',
            'my.table$%s' % name,
            xlab = '"%s"' % name,
            ylab = '"counts"',
            main='"%s for %d selection samples"' % (
                name, nsels))
    user_script_content = out.getvalue()
    #j
    #print 'table string:'
    #print table_string
    #print 'user script:'
    #print user_script_content
    #
    # call R to get the tikz code
    retcode, r_out, r_err, tikz_code = RUtil.run_plotter(
            table_string, user_script_content, 'tikz',
            width=5, height=5)
    if retcode:
        raise RUtil.RError(r_err)
    return tikz_code

def get_latex_documentbody(fs):
    """
    The latex documentbody should have a bunch of tikz pieces in it.
    Each tikz piece should have been generated from R.
    """
    nstates = fs.nresidues ** fs.nsites
    if nstates > 256:
        raise ValueError('the mutation rate matrix is too big')
    # get the mutation matrix
    if fs.sparse:
        Q_mut = mrate.get_sparse_sequence_rate_matrix(fs.nresidues, fs.nsites)
    elif fs.dense:
        Q_mut = mrate.get_dense_sequence_rate_matrix(fs.nresidues, fs.nsites)
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
    # compute the statistics
    ER_ratios, NSR_ratios, ER_NSR_ratios  = get_statistic_ratios(Q_mut, Q_sels)
    M = zip(*(ER_ratios, NSR_ratios, ER_NSR_ratios))
    column_headers = ('ER.ratio', 'NSR.ratio', 'ER.times.NSR.ratio')
    table_string = RUtil.get_table_string(M, column_headers)
    nsels = len(Q_sels)
    # write the latex code
    out = StringIO()
    for name in column_headers:
        print >> out, get_r_tikz_hist(nsels, table_string, name)
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
