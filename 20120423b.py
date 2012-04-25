"""
Plot log(ratio(E(L))) and E(log(ratio(L))) for non-uniform stationary process.
"""

from StringIO import StringIO
import math

import numpy as np

import Form
import FormOut
import mrate
import ctmcmi
import MatrixUtil
import RUtil
from RUtil import mk_call_str

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.Sequence('lowtri',
                'strictly lower triangular mut exch',
                ('1', '1 1', '1 1 1')),
            Form.Sequence('mutweights',
                'unnormalized mut stationary distn',
                ('1', '1', '1', '1')),
            Form.Float('mutscale',
                'extra mutation process scaling factor',
                '1.333333333333', low_exclusive=0),
            Form.Sequence('mutselweights',
                'unnormalized mut-sel stationary distn',
                ('1', '1', '0.001', '0.001')),
            Form.ImageFormat()]
    return form_objects

def get_form_out():
    return FormOut.Image('mutual-info-and-analog')

def get_presets():
    presets = [
            Form.Preset(
                'only mutual info crosses',
                {
                    'lowtri' : ('1', '0 1', '1 0 1'),
                    'mutweights' : ('1', '0.001', '1', '0.001'),
                    'mutscale' : '8',
                    'mutselweights' : ('0.7', '0.1', '0.1', '0.1'),
                    'imageformat' : 'png'}),
            Form.Preset(
                'only the analog crosses',
                {
                    'lowtri' : ('1', '1 1', '1 1 1'),
                    'mutweights' : ('1', '0.001', '1', '0.001'),
                    'mutscale' : '0.5',
                    'mutselweights' : ('0.7', '0.1', '0.1', '0.1'),
                    'imageformat' : 'png'})]
    return presets

def get_response_content(fs):
    M, R = get_input_matrices(fs)
    # create the R table string and scripts
    headers = [
            't',
            'mi.true.mut',
            'mi.true.mutsel',
            'mi.analog.mut',
            'mi.analog.mutsel']
    npoints = 100
    t_low = 0.0
    t_high = 5.0
    t_incr = (t_high - t_low) / (npoints - 1)
    t_values = [t_low + t_incr*i for i in range(npoints)]
    # get the data for the R table
    arr = []
    for t in t_values:
        mi_mut = ctmcmi.get_mutual_information(M, t)
        mi_mutsel = ctmcmi.get_mutual_information(R, t)
        mi_analog_mut = ctmcmi.get_ll_ratio_wrong(M, t)
        mi_analog_mutsel = ctmcmi.get_ll_ratio_wrong(R, t)
        row = [t, mi_mut, mi_mutsel, mi_analog_mut, mi_analog_mutsel]
        arr.append(row)
    # get the R table
    table_string = RUtil.get_table_string(arr, headers)
    # get the R script
    script = get_ggplot()
    # create the R plot image
    device_name = Form.g_imageformat_to_r_function[fs.imageformat]
    retcode, r_out, r_err, image_data = RUtil.run_plotter(
            table_string, script, device_name)
    if retcode:
        raise RUtil.RError(r_err)
    return image_data

def get_plainplot():
    out = StringIO()
    print >> out, mk_call_str(
            'plot',
            'my.table$t',
            'my.table$mutual.info.mutsel',
            type='"n"',
            xlab='"time"',
            ylab='""',
            ylim=mk_call_str('c', '0', '1.0'),
            main='"MI (black) and MI analog (red) over time"')
    print >> out, mk_call_str(
            'lines',
            'my.table$t',
            'my.table$mutual.info.mut',
            col='"black"')
    print >> out, mk_call_str(
            'lines',
            'my.table$t',
            'my.table$mutual.info.analog.mut',
            col='"red"')
    print >> out, mk_call_str(
            'lines',
            'my.table$t',
            'my.table$mutual.info.mutsel',
            col='"black"')
    print >> out, mk_call_str(
            'lines',
            'my.table$t',
            'my.table$mutual.info.analog.mutsel',
            col='"red"')
    return out.getvalue()

def get_ggplot():
    out = StringIO()
    print >> out, mk_call_str('require', '"reshape"')
    print >> out, mk_call_str('require', '"ggplot2"')
    print >> out, 'my.table.long <-',
    print >> out, mk_call_str('melt', 'my.table', id='"t"')
    print >> out, 'ggplot(data=my.table.long,'
    print >> out, mk_call_str('aes', x='t', y='value', colour='variable')
    print >> out, ') + geom_line()',
    print >> out, '+',
    print >> out, mk_call_str(
            'xlim', '0',
            mk_call_str('max', 'my.table.long$t')),
    print >> out, '+',
    print >> out, mk_call_str(
            'ylim', '0',
            mk_call_str('max', 'my.table.long$value'))
    return out.getvalue()

def get_input_matrices(fs):
    """
    @return: M, R
    """
    # get the positive strict lower triangular part of the S matrix
    L = []
    for i, line in enumerate(fs.lowtri):
        values = line.split()
        if len(values) != i + 1:
            msg = 'expected %d values on line "%s"' % (
                    i+1, line)
            raise ValueError(msg)
        vs = [float(v) for v in values]
        if any(x<0 for x in vs):
            raise ValueError('exchangeabilities must be nonnegative')
        L.append(vs)
    # get the mut and mutsel weights
    mut_weights = [float(v) for v in fs.mutweights]
    mutsel_weights = [float(v) for v in fs.mutselweights]
    if any(x<=0 for x in mut_weights + mutsel_weights):
        raise ValueError('stationary weights must be positive')
    # normalize weights to distributions
    mut_distn = [v / sum(mut_weights) for v in mut_weights]
    mutsel_distn = [v / sum(mutsel_weights) for v in mutsel_weights]
    # get the exchangeability matrix
    nstates = len(L) + 1
    S = np.zeros((nstates, nstates))
    for i, row in enumerate(L):
        for j, v in enumerate(row):
            S[i+1, j] = v
            S[j, i+1] = v
    # check the state space sizes implied by the inputs
    if len(set(len(x) for x in (S, mut_weights, mutsel_weights))) != 1:
        msg = 'the inputs do not agree on the state space size'
        raise ValueError(msg)
    # check for sufficient number of states
    if nstates < 2:
        msg = 'at least two states are required'
        raise ValueError(msg)
    # check reducibility of the exchangeability
    MatrixUtil.assert_symmetric_irreducible(S)
    # get the mutation rate matrix
    M = S * mut_distn * fs.mutscale
    M -= np.diag(np.sum(M, axis=1))
    # check sign symmetry and irreducibility
    MatrixUtil.assert_symmetric_irreducible(np.sign(M))
    # get the mutation selection balance rate matrix
    R = mrate.to_gtr_halpern_bruno(M, mutsel_distn)
    # check sign symmetry and irreducibility
    MatrixUtil.assert_symmetric_irreducible(np.sign(R))
    # check the stationary distributions
    mut_distn_observed = mrate.R_to_distn(M)
    if not np.allclose(mut_distn_observed, mut_distn):
        msg = 'internal mut stationary distribution computation error'
        raise ValueError(msg)
    mutsel_distn_observed = mrate.R_to_distn(R)
    if not np.allclose(mutsel_distn_observed, mutsel_distn):
        msg = 'internal mut-sel stationary distribution computation error'
        raise ValueError(msg)
    # return the values
    return M, R

