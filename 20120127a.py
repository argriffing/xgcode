"""
Write explicit bounds for CTMC mutual information asymptotics.

By CTMC I mean
a reversible irreducible finite-state continuous-time Markov chain.
By bounds for mutual information asymptotics I mean
upper and lower bounds on the mutual information between two
observations separated by time \\( t \\),
when \\( t \\) is large enough.
The specific upper and lower bounds are each
positive decreasing pure exponential functions that decay towards zero
at the same decay rate, but with different coefficients.
The decay rate of the bounding exponential functions
is an eigenvalue of the rate matrix.
"""


from StringIO import StringIO
import math
import random

import numpy as np
import scipy
from scipy import optimize

import Form
import FormOut
import mrate
import ctmcmi
import ctmcmitaylor



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
            Form.RadioGroup('sel_var', 'selection parameter variance', [
                Form.RadioItem('really_low_var', 'really low variance'),
                Form.RadioItem('low_var', 'low variance'),
                Form.RadioItem('medium_var', 'medium variance', True),
                Form.RadioItem('high_var', 'high variance'),
                Form.RadioItem('really_high_var', 'really high variance')]),
            Form.RadioGroup('sel_skew', 'selection parameter skew', [
                Form.RadioItem('neg_skew', 'some are very fit'),
                Form.RadioItem('no_skew', 'no skew', True),
                Form.RadioItem('pos_skew', 'some are very unfit')])]
    return form_objects

def get_form_out():
    return FormOut.Report()

def sample_rate_matrix(fs, Q_mut):
    nstates = len(Q_mut)
    # sample the selection parameters
    if fs.really_low_var:
        v = 0.04
    elif fs.low_var:
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
    return Q

class RateProperties:
    def __init__(self, Q):
        self.Q = Q
        self.relaxation_time = mrate.R_to_relaxation_time(Q)
        self.p = min(mrate.R_to_distn(Q))
        self.N = len(Q)
        self.lam = - 1 / self.relaxation_time
        key_time_points = ctmcmitaylor.get_key_time_points(
            self.lam, self.p, self.N)
        self.time_to_uniformity, self.time_to_usefulness = key_time_points
    def __str__(self):
        out = StringIO()
        print >> out, 'rate matrix:'
        print >> out, self.Q
        print >> out
        print >> out, 'relaxation time:'
        print >> out, self.relaxation_time
        print >> out
        print >> out, 'min stationary probability:'
        print >> out, self.p
        print >> out
        print >> out, 'expected rate:'
        print >> out, mrate.R_to_total_rate(self.Q)
        print >> out
        print >> out, 'time to uniform-over-entries taylor h_2 bound:'
        print >> out, self.time_to_uniformity
        print >> out
        print >> out, 'time to informativeness for weak inequality:'
        print >> out, self.time_to_usefulness
        print >> out
        return out.getvalue().rstrip()

def get_response_content(fs):
    nstates = fs.nresidues ** fs.nsites
    if nstates > 256:
        raise ValueError('the mutation rate matrix is too big')
    # get the mutation matrix
    Q_mut = mrate.get_sparse_sequence_rate_matrix(fs.nresidues, fs.nsites)
    # get the random selection matrix which we will use from now on
    Q_sel = sample_rate_matrix(fs, Q_mut)
    # define the time points
    #incr = (fs.t_high - fs.t_low) / (fs.ntimes - 1)
    #times = [fs.t_low + i*incr for i in range(fs.ntimes)]
    mut_info = RateProperties(Q_mut)
    sel_info = RateProperties(Q_sel)
    # compute the intersection time
    x_time_top = math.log(2 * nstates - 1)
    x_time_bot = 2 * abs(mut_info.lam - sel_info.lam)
    x_time = x_time_top / x_time_bot
    # compute the upper bound on the judgement time
    T_second_order = max(
            x_time,
            mut_info.time_to_usefulness,
            sel_info.time_to_usefulness)
    # define the name of the eventually winning process
    if mut_info.relaxation_time > sel_info.relaxation_time:
        x = 'mutation'
        slow_info = mut_info
        fast_info = sel_info
    else:
        x = 'mutation-selection balance'
        slow_info = sel_info
        fast_info = mut_info
    eventual_winner_name = 'the %s process' % x
    # get a more sophisticated bound
    third_order_x_time = ctmcmitaylor.get_sophisticated_time_bound(
            -slow_info.lam,
            -fast_info.lam,
            slow_info.N,
            fast_info.N,
            slow_info.p,
            fast_info.p)
    if third_order_x_time is not None:
        T_third_order = max(
                third_order_x_time,
                mut_info.time_to_uniformity,
                sel_info.time_to_uniformity)
    else:
        T_third_order = None
    # Define a naive crossing time.
    # This is not a bound on the true mutual information doomsday,
    # but it shows a limit of our approach.
    # It is the bound on the spectral taylor approximation
    # given only the second eigenvalues and not the other ones.
    naive_x_time_top = math.log(nstates - 1)
    naive_x_time_bot = 2 * abs(mut_info.lam - sel_info.lam)
    naive_x_time = naive_x_time_top / naive_x_time_bot
    # write the report
    np.set_printoptions(linewidth=200)
    out = StringIO()
    print >> out, '*** mutation rate matrix info ***'
    print >> out
    print >> out, mut_info
    print >> out
    print >> out
    print >> out, '*** mutation-selection balance rate matrix info ***'
    print >> out
    print >> out, sel_info
    print >> out
    print >> out
    print >> out, '*** note ***'
    print >> out
    print >> out, 'with the general approach taken here,'
    print >> out, 'we will not find an eigenvalue time bound'
    print >> out, 'smaller than', naive_x_time
    print >> out
    print >> out
    print >> out, '*** weak inequality ***'
    print >> out
    print >> out, 'When t >', T_second_order, eventual_winner_name
    print >> out, 'has greater mutual information (MI) and approximate MI.'
    print >> out
    print >> out
    print >> out, '*** stronger inequality ***'
    print >> out
    if T_third_order is None:
        print >> out, 'the numerical solver failed to converge'
    else:
        print >> out, 'When t >', T_third_order, eventual_winner_name
        print >> out, 'has greater mutual information (MI) and approximate MI.'
    return out.getvalue()

