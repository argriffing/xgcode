r"""
Plot max information over all selection parameters.

For an interval of time, the divergence time information is plotted
for each of three processes.
The first process is site-independent neutral evolution.
The second process is the site-independent mutation-selection balance
process whose selection parameters
give the greatest information at the given time.
The third process is the site-dependent mutation-selection balance
process whose selection parameters
give the greatest information at the given time.
Many of the functions associated with this web script
have been refitted to accept a stationary distribution hint
instead of having to compute the stationary distribution themselves
from the rate matrix.
Presumably this helps stability as the processes approach reducibility.
Maybe this could be misleading if there is a discontinuity
as the process becomes reducible.
For the purposes of this plot, site dependence means the non-independent
evolution of two nucleotide positions in a DNA alignment.
"""

from StringIO import StringIO
import argparse
import math
import time
import random
import heapq
from itertools import product

import numpy as np
import scipy
from scipy import linalg
from scipy import optimize

import Form
import FormOut
import ctmcmi
import mrate
import divtime
import cheeger
import MatrixUtil
from MatrixUtil import ndot
import RUtil
from RUtil import mk_call_str

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.Sequence('lowtri',
                'stricty lower triangular mut exch',
                ('1',)),
            Form.Sequence('distn_weights',
                'unnormalized mut stationary distn',
                ('1', '1')),
            Form.Float('scale', 'extra scaling factor',
                '1.0', low_exclusive=0),
            Form.Float('start_time', 'start time', '0.1', low_exclusive=0),
            Form.Float('stop_time', 'stop time', '0.6', low_exclusive=0),
            Form.Integer('nsites', 'number of sites',
                '2', low=2, high=4),
            Form.CheckGroup('proctypes', 'foo types', [
                Form.CheckItem('indep_mutation',
                    'site independent mutation', True),
                Form.CheckItem('indep_balance',
                    'site independent mutation-selection balance', True),
                Form.CheckItem('dep_balance',
                    'site dependent mutation-selection balance', True)]),
            Form.RadioGroup('infotype', 'information variant', [
                Form.RadioItem('info_fis', 'Fisher information'),
                Form.RadioItem('info_mut', 'mutual information', True)]),
            Form.ImageFormat(),
            ]
    return form_objects

def get_form_out():
    return FormOut.Image('plot')

def gen_overdispersed_events(low, high):
    """
    This is a generator that samples overdispersed events in an interval.
    It is useful for plotting.
    """
    # Create the min heap.
    # The triples are the neg length, the low, and the high values.
    # Popping from the queue will return the longest interval.
    # The queue grows linearly with the number of sampled events.
    yield low
    yield high
    q = [(low-high, low, high)]
    while True:
        dummy, a, b = heapq.heappop(q)
        mid = random.uniform(a, b)
        heapq.heappush(q, (a-mid, a, mid))
        heapq.heappush(q, (mid-b, mid, b))
        yield mid

def get_site_independent_distn(v, nsites):
    n = len(v)
    v_out = np.zeros(n**nsites)
    for values in product(range(n), repeat=nsites):
        index = sum(k * (n**i) for i, k in enumerate(values))
        p = 1
        for value in values:
            p *= v[value]
        v_out[index] = p
    return v_out

def get_site_independent_process(R, nsites):
    """
    @param R: single site process
    @param nsites: return a rate matrix that acts on this many sites
    @return: a site independent rate matrix
    """
    n = len(R)
    nprod = n**nsites
    Q = np.zeros((nprod, nprod))
    for a_values in product(range(n), repeat=nsites):
        a_index = sum(k * (n**i) for i, k in enumerate(a_values))
        for b_values in product(range(n), repeat=nsites):
            b_index = sum(k * (n**i) for i, k in enumerate(b_values))
            diffs = [(i, a, b) for i, (a, b) in enumerate(zip(
                a_values, b_values)) if a != b]
            if len(diffs) == 1:
                i, a, b = diffs[0]
                Q[a_index, b_index] = R[a, b]
    Q -= np.diag(np.sum(Q, axis=1))
    return Q

class OptDep:
    """
    This is for numerical optimization.
    """
    def __init__(self, R, v, t, f_info, f_selection):
        """
        @param R: two site independent mutation rate matrix
        @param v: stationary distribution of R
        @param t: time
        @param f_info: info function that takes a rate matrix and a time
        @param f_selection: takes a rate matrix and a current and target distn
        """
        self.R = R
        self.v = v
        self.t = t
        self.f_info = f_info
        self.f_selection = f_selection
    def get_process(self, X):
        """
        @param X: some log ratio probabilities
        @return: the mut-sel balance rate matrix and stationary distn
        """
        if len(X) != len(self.R) - 1:
            raise ValueError('state space size mismatch')
        # define the mutation and balance energies
        u = -(np.log(self.v) - math.log(self.v[0]))
        v = np.array([0.0] + X.tolist())
        # apply the selection using the stable halpern bruno function
        Q = self.f_selection(self.R, u, v)
        # get the stationary distribution for downstream
        qdistn = np.exp(-v)
        qdistn /= np.sum(qdistn)
        if np.any(np.isnan(Q)):
            print self.R
            print qdistn
            print Q
            raise ValueError('the rate matrix has nans')
        return Q, qdistn
    def __call__(self, X):
        """
        @param X: some log ratio probabilities
        @return: neg info value for minimization
        """
        Q, v_target = self.get_process(X)
        return -self.f_info(Q, v_target, self.t)

class OptIndep:
    """
    This is for numerical optimization.
    """
    def __init__(self, R, v, nsites, t, f_info, f_selection):
        """
        @param R: single state mutation rate matrix
        @param v: stationary distribution of R
        @parma nsites: number of independently evolving sites
        @param t: time
        @param f_info: info function that takes a rate matrix and a time
        @param f_selection: takes a rate matrix and a current and target distn
        """
        self.R = R
        self.v = v
        self.nsites = nsites
        self.t = t
        self.f_info = f_info
        self.f_selection = f_selection
    def get_process(self, X):
        """
        @param X: some energies defining the probability distribution
        @return: the mut-sel balance rate matrix and stationary distn
        """
        if len(X) != len(self.R) - 1:
            raise ValueError('state space size mismatch')
        # define the mutation and balance energies
        u = -(np.log(self.v) - math.log(self.v[0]))
        v = np.array([0.0] + X.tolist())
        # apply the selection using the stable halpern bruno function
        Q_single_site = self.f_selection(self.R, u, v)
        Q = get_site_independent_process(Q_single_site, self.nsites)
        # get the stationary distribution for downstream
        qdistn = np.exp(-v)
        qdistn /= np.sum(qdistn)
        qdistn_site_indep = get_site_independent_distn(qdistn, self.nsites)
        if np.any(np.isnan(Q)):
            print self.R
            print qdistn_site_indep
            print Q
            raise ValueError('the rate matrix has nans')
        return Q, qdistn_site_indep
    def __call__(self, X):
        """
        @param X: some log ratio probabilities
        @return: neg info value for minimization
        """
        Q, v_site_indep = self.get_process(X)
        return -self.f_info(Q, v_site_indep, self.t)

def get_response_content(fs):
    # hardcode the amount of time allowed
    nseconds = 4
    # validate and store user input
    if fs.stop_time <= fs.start_time:
        raise ValueError('check the start and stop times')
    M = get_input_matrix(fs)
    nstates = len(M)
    nsites = fs.nsites
    if nstates ** nsites > 16:
        raise ValueError('the site dependent rate matrix is too big')
    # precompute some stuff
    M_site_indep = get_site_independent_process(M, nsites)
    v = mrate.R_to_distn(M)
    v_site_indep = get_site_independent_distn(v, nsites)
    if fs.info_fis:
        f_info = divtime.get_fisher_info_known_distn_fast
    elif fs.info_mut:
        f_info = ctmcmi.get_mutual_info_known_distn
    else:
        raise ValueError('no info type specified')
    f_selection = mrate.to_gtr_hb_known_energies
    # Spend a lot of time doing the optimizations
    # to construct the points for the R table.
    t0 = time.time()
    arr = []
    for t in gen_overdispersed_events(fs.start_time, fs.stop_time):
        if time.time() - t0 > nseconds:
            break
        row = [t]
        # get the site-dependent mutation selection balance information
        if fs.dep_balance:
            dep_balance = OptDep(
                    M_site_indep, v_site_indep, t, f_info, f_selection)
            X0 = np.random.randn(nstates ** nsites - 1)
            xopt = scipy.optimize.fmin(dep_balance, X0)
            max_dep_balance_info = -dep_balance(xopt)
            row.append(max_dep_balance_info)
            # for debug
            Q_bal, v_bal = dep_balance.get_process(xopt)
            print 'dependent balance:'
            print max_dep_balance_info
            print v_bal
            print Q_bal
        # get the site-independent mutation selection balance information
        if fs.indep_balance:
            indep_balance = OptIndep(
                    M, v, nsites, t, f_info, f_selection)
            X0 = np.random.randn(nstates-1)
            xopt = scipy.optimize.fmin(indep_balance, X0)
            max_indep_balance_info = -indep_balance(xopt)
            row.append(max_indep_balance_info)
            # for debug
            Q_bal, v_bal = indep_balance.get_process(xopt)
            print 'independent balance:'
            print max_indep_balance_info
            print v_bal
            print Q_bal
        # get the site-independent mutation process information
        if fs.indep_mutation:
            indep_mut_info = f_info(M_site_indep, v_site_indep, t)
            row.append(indep_mut_info)
        # add the data row to the table
        arr.append(row)
    arr.sort()
    npoints = len(arr)
    # create the R table string and scripts
    headers = ['t']
    if fs.dep_balance:
        headers.append('max.site.dep.balance')
    if fs.indep_balance:
        headers.append('max.site.indep.balance')
    if fs.indep_mutation:
        headers.append('site.indep.mutation')
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
            msg = 'expected %d values on line "%s"' % (
                    i+1, line)
            raise ValueError(msg)
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
        msg = 'the inputs do not agree on the state space size'
        raise ValueError(msg)
    # check for sufficient number of states
    if nstates < 2:
        msg = 'at least two states are required'
        raise ValueError(msg)
    # check reducibility of the exchangeability
    if not MatrixUtil.is_symmetric_irreducible(S):
        raise ValueError('exchangeability is not irreducible')
    # get the mutation rate matrix
    M = S * distn * fs.scale
    M -= np.diag(np.sum(M, axis=1))
    # check sign symmetry and irreducibility
    if not MatrixUtil.is_symmetric_irreducible(np.sign(M)):
        msg = 'mutation rate matrix is not sign symmetric irreducible'
        raise ValueError(msg)
    # check the stationary distributions
    distn_observed = mrate.R_to_distn(M)
    if not np.allclose(distn_observed, distn):
        msg = 'internal mut stationary distribution computation error'
        raise ValueError(msg)
    # return the values
    return M

