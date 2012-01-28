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
    # mutual information proportion
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
    v4 = [math.exp(-t*mrate.R_to_total_rate(Q)) for Q in Q_sels]
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
    # add the mutual information sign proportion
    statistics.append(sum(1 for x in mi_signs if x == 1) / float(nsels))
    return statistics

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


def get_time_stats(fs):
    nselections = 1
    nstates = fs.nresidues ** fs.nsites
    if nstates > 256:
        raise ValueError('the mutation rate matrix is too big')
    # get the mutation matrix
    Q_mut = mrate.get_sparse_sequence_rate_matrix(fs.nresidues, fs.nsites)
    # sample a bunch of mutation-selection rate matrices
    Q_sels = [sample_rate_matrix(fs, Q_mut) for i in range(nselections)]
    # define the time points
    incr = (fs.t_high - fs.t_low) / (fs.ntimes - 1)
    times = [fs.t_low + i*incr for i in range(fs.ntimes)]
    # compute the statistics
    time_stats = [get_time_point_summary(Q_mut, Q_sels, t) for t in times]
    return time_stats

def get_script_a(divname, transposed_time_stats):
    """
    A summary of mutation-selection mutual information.
    @param transposed_time_stats: statistics on a timeline
    @return: a javascript that uses dygraph to display the timeline
    """
    # define the containing javascript div
    containing_div = 'document.getElementById("%s")' % divname
    # define the csv as a javascript string
    header = 'time,mut,sel'
    ts, mis, mins, lows, means, highs, maxs = transposed_time_stats[:7]
    pairs = ['%0.4f,%s;%s;%s,%s;%s;%s' % v for v in zip(*(
        ts,mis,mis,mis,mins,means,maxs))]
    csv_string = '"' + '\\n'.join([header] + pairs) + '\\n"'
    option_string = '{customBars: true}'
    # Write the script; it has only one call.
    return 'g = new Dygraph(%s, %s, %s);' % (
            containing_div, csv_string, option_string)

def get_script_b(divname, transposed_time_stats):
    """
    Probability that selection mutual information is bigger.
    @param transposed_time_stats: statistics on a timeline
    @return: a javascript that uses dygraph to display the timeline
    """
    # define the containing javascript div
    containing_div = 'document.getElementById("%s")' % divname
    # define the csv as a javascript string
    header = 'time,prop'
    ts = transposed_time_stats[0]
    props = transposed_time_stats[-1]
    pairs = ['%0.4f,%s' % v for v in zip(*(ts, props))]
    csv_string = '"' + '\\n'.join([header] + pairs) + '\\n"'
    option_string = '{}'
    # Write the script; it has only one call.
    return 'g = new Dygraph(%s, %s, %s);' % (
            containing_div, csv_string, option_string)

def get_script_c(divname, transposed_time_stats):
    """
    Correlation with MI.
    @param transposed_time_stats: statistics on a timeline
    @return: a javascript that uses dygraph to display the timeline
    """
    # define the containing javascript div
    containing_div = 'document.getElementById("%s")' % divname
    # define the csv as a javascript string
    header = 'time,fn1,fn2,fn3,fn4,fn5'
    ts = transposed_time_stats[0]
    corr1, corr2, corr3, corr4, corr5 = transposed_time_stats[7:7+5]
    props = transposed_time_stats[-1]
    pairs = ['%0.4f,%s,%s,%s,%s,%s' % v for v in zip(*(
        ts, corr1, corr2, corr3, corr4, corr5))]
    csv_string = '"' + '\\n'.join([header] + pairs) + '\\n"'
    option_string = '{}'
    # Write the script; it has only one call.
    return 'g = new Dygraph(%s, %s, %s);' % (
            containing_div, csv_string, option_string)

def get_script_d(divname, transposed_time_stats):
    """
    Proportion same sign.
    @param transposed_time_stats: statistics on a timeline
    @return: a javascript that uses dygraph to display the timeline
    """
    # define the containing javascript div
    containing_div = 'document.getElementById("%s")' % divname
    # define the csv as a javascript string
    header = 'time,fn1,fn2,fn3,fn4,fn5'
    ts = transposed_time_stats[0]
    prop1, prop2, prop3, prop4, prop5 = transposed_time_stats[12:12+5]
    props = transposed_time_stats[-1]
    pairs = ['%0.4f,%s,%s,%s,%s,%s' % v for v in zip(*(
        ts, prop1, prop2, prop3, prop4, prop5))]
    csv_string = '"' + '\\n'.join([header] + pairs) + '\\n"'
    option_string = '{}'
    # Write the script; it has only one call.
    return 'g = new Dygraph(%s, %s, %s);' % (
            containing_div, csv_string, option_string)

def get_script_e(divname, transposed_time_stats):
    """
    Info wrt mutual info.
    @param transposed_time_stats: statistics on a timeline
    @return: a javascript that uses dygraph to display the timeline
    """
    # define the containing javascript div
    containing_div = 'document.getElementById("%s")' % divname
    # define the csv as a javascript string
    header = 'time,fn1,fn2,fn3,fn4,fn5'
    ts = transposed_time_stats[0]
    info1, info2, info3, info4, info5 = transposed_time_stats[17:17+5]
    props = transposed_time_stats[-1]
    pairs = ['%0.4f,%s,%s,%s,%s,%s' % v for v in zip(*(
        ts, info1, info2, info3, info4, info5))]
    csv_string = '"' + '\\n'.join([header] + pairs) + '\\n"'
    option_string = '{}'
    # Write the script; it has only one call.
    return 'g = new Dygraph(%s, %s, %s);' % (
            containing_div, csv_string, option_string)

def get_dygraph_html(script_fn, divname, time_stats):
    """
    @param divname: the element id of the containing div
    @param script_fn: get the script given the element id and time stats
    @return: a chunk of html
    """
    out = StringIO()
    print >> out, '<div id="%s"></div>' % divname
    print >> out, '<script type="text/javascript">'
    print >> out, script_fn(divname, time_stats)
    print >> out, '</script>'
    return out.getvalue().rstrip()

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
        print >> out, 'time to uniformity:'
        print >> out, self.time_to_uniformity
        print >> out
        print >> out, 'time to usefulness:'
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
    T = max(x_time, mut_info.time_to_usefulness, sel_info.time_to_usefulness)
    # define the name of the eventually winning process
    if mut_info.relaxation_time > sel_info.relaxation_time:
        x = 'mutation'
    else:
        x = 'mutation-selection balance'
    eventual_winner_name = 'the %s process' % x
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
    print >> out, '*** asymptotic inequality ***'
    print >> out
    print >> out, 'When t > %s %s has greater mutual information.' % (
            T, eventual_winner_name)
    return out.getvalue()

