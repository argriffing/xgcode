"""
Plot upper and lower mutual information bounds for reversible Markov processes.

The nonspectral upper bound is not so great; it is just the constant entropy
of the stationary distribution.
There are two nonspectral lower bounds, each of which corresponds to the
mutual information of a process which I believe never has more
mutual information than the original process.
One of these less informative processes
is an F81-ization which has the same stationary distribution so it
is a tight lower bound at times near zero.
The other less informative process
is a two state process defined by the strongest barrier
in the original process;
it is a weak lower bound at small divergence times
because its stationary distribution has not so much entropy,
but at larger times it should be a quite good lower bound.
I have not proved that the spectral bounds are actually bounds,
but they are close approximations near zero and infinity.
For that matter I have not proved that the non-spectral
bounds are bounds either.
"""

from StringIO import StringIO
import math

import numpy as np

import Form
import FormOut
import mrate
import ctmcmi
import msimpl
import MatrixUtil
import RUtil
from RUtil import mk_call_str

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.Sequence('lowtri',
                'strictly lower triangular exchangeabilities',
                ('1', '1 1', '1 0 0')),
            Form.Sequence('distn_weights',
                'unnormalized stationary distribution',
                ('1', '0.1', '0.1', '1')),
            Form.Float('scale',
                'extra scaling factor',
                '1.333333333333', low_exclusive=0),
            Form.FloatInterval(
                'start_time', 'stop_time', 'divtime interval',
                '0', '5', low_inclusive=0, low_width_exclusive=0),
            Form.CheckGroup('options', 'options', [
                Form.CheckItem('show_entropy', 'show entropy', True)]),
            Form.ImageFormat()]
    return form_objects

def get_form_out():
    return FormOut.Image('mutual-info-bounds')

def get_response_content(fs):
    M = get_input_matrix(fs)
    # create the R table string and scripts
    headers = ['t']
    if fs.show_entropy:
        headers.append('ub.entropy')
    headers.extend([
            'ub.jc.spectral',
            'ub.f81.spectral',
            'mutual.information',
            'lb.2.state.spectral',
            'lb.2.state',
            'lb.f81',
            ])
    npoints = 100
    t_low = fs.start_time
    t_high = fs.stop_time
    t_incr = (t_high - t_low) / (npoints - 1)
    t_values = [t_low + t_incr*i for i in range(npoints)]
    # define some extra stuff
    v = mrate.R_to_distn(M)
    entropy = -np.dot(v, np.log(v))
    n = len(M)
    gap = sorted(abs(x) for x in np.linalg.eigvals(M))[1]
    print 'stationary distn:', v
    print 'entropy:', entropy
    print 'spectral gap:', gap
    M_slow_jc = gap * (1.0 / n) * (np.ones((n,n)) - n*np.eye(n))
    M_slow_f81 = gap * np.outer(np.ones(n), v)
    M_slow_f81 -= np.diag(np.sum(M_slow_f81, axis=1))
    M_f81 = msimpl.get_fast_f81(M)
    M_2state = msimpl.get_fast_two_state_autobarrier(M)
    M_2state_spectral = -gap * M_2state / np.trace(M_2state)
    # get the data for the R table
    arr = []
    for u in t_values:
        # experiment with log time
        #t = math.exp(u)
        t = u
        mi_slow_jc = ctmcmi.get_mutual_information(M_slow_jc, t)
        mi_slow_f81 = ctmcmi.get_mutual_information(M_slow_f81, t)
        mi_mut = ctmcmi.get_mutual_information(M, t)
        mi_2state_spectral = ctmcmi.get_mutual_information(M_2state_spectral, t)
        mi_f81 = ctmcmi.get_mutual_information(M_f81, t)
        mi_2state = ctmcmi.get_mutual_information(M_2state, t)
        row = [u]
        if fs.show_entropy:
            row.append(entropy)
        row.extend([mi_slow_jc, mi_slow_f81,
                mi_mut, mi_2state_spectral, mi_2state, mi_f81])
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
            'xlim',
            mk_call_str('min', 'my.table.long$t'),
            mk_call_str('max', 'my.table.long$t')),
    print >> out, '+',
    print >> out, mk_call_str(
            'ylim', '0',
            mk_call_str('max', 'my.table.long$value'))
    return out.getvalue()

def get_input_matrix(fs):
    """
    @return: M
    """
    # get the positive strict lower triangular part of the S matrix
    L = []
    for i, line in enumerate(fs.lowtri):
        values = line.split()
        if len(values) != i + 1:
            raise ValueError(
                    'expected %d values on line "%s"' % (i+1, line))
        vs = [float(v) for v in values]
        if any(x<0 for x in vs):
            raise ValueError('exchangeabilities must be nonnegative')
        L.append(vs)
    # get the stationary distribution weights
    distn_weights = [float(v) for v in fs.distn_weights]
    if any(x<=0 for x in distn_weights):
        raise ValueError('stationary weights must be positive')
    # normalize weights to distributions
    distn = [v / sum(distn_weights) for v in distn_weights]
    # get the exchangeability matrix
    nstates = len(L) + 1
    S = np.zeros((nstates, nstates))
    for i, row in enumerate(L):
        for j, v in enumerate(row):
            S[i+1, j] = v
            S[j, i+1] = v
    # check the state space sizes implied by the inputs
    if len(set(len(x) for x in (S, distn_weights))) != 1:
        raise ValueError('the inputs do not agree on the state space size')
    # check for sufficient number of states
    if nstates < 2:
        raise ValueError('at least two states are required')
    # check reducibility of the exchangeability
    if not MatrixUtil.is_symmetric_irreducible(S):
        raise ValueError('exchangeability is not irreducible')
    # get the mutation rate matrix
    M = S * distn * fs.scale
    M -= np.diag(np.sum(M, axis=1))
    # check sign symmetry and irreducibility
    if not MatrixUtil.is_symmetric_irreducible(np.sign(M)):
        raise ValueError(
                'mutation rate matrix is not sign symmetric irreducible')
    # check the stationary distributions
    distn_observed = mrate.R_to_distn(M)
    if not np.allclose(distn_observed, distn):
        raise ValueError(
                'internal mut stationary distribution computation error')
    # return the values
    return M

