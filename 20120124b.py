"""
Plot mutual information using the dygraph javascript.
"""


from StringIO import StringIO
import math
import random

import numpy as np

import Form
import FormOut
import mrate
import ctmcmi


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
            Form.Float('t_low', 'initial time',
                '0.001', low_exclusive=0),
            Form.Float('t_high', 'final time',
                '4.0', low_exclusive=0),
            Form.Integer('ntimes', 'sample this many time points',
                10, low=3, high=100)]
    return form_objects

def get_form_out():
    return FormOut.Html('timeline')

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
    return statistics

def get_time_stats(fs):
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
    time_stats = [get_time_point_summary(Q_mut, Q_sels, t) for t in times]
    return time_stats

def get_script(time_stats):
    """
    @param time_stats: statistics on a timeline
    @return: a javascript that uses dygraph to display the timeline
    """
    # define the containing javascript div
    containing_div = 'document.getElementById("graphdiv")'
    # define the csv as a javascript string
    header = 'time,mut,sel'
    ts, mis, mins, lows, means, highs, maxs = zip(*time_stats)
    pairs = ['%0.4f,%s;%s;%s,%s;%s;%s' % v for v in zip(*(
        ts,mis,mis,mis,mins,means,maxs))]
    csv_string = '"' + '\\n'.join([header] + pairs) + '\\n"'
    option_string = '{customBars: true}'
    # Write the script; it has only one call.
    return 'g = new Dygraph(%s, %s, %s);' % (
            containing_div, csv_string, option_string)

def get_response_content(fs):
    dygraphs_url = 'http://dygraphs.com/dygraph-combined.js'
    time_stats = get_time_stats(fs)
    script = get_script(time_stats)
    # create an html file that uses the javascript link
    out = StringIO()
    print >> out, '<html>'
    print >> out, '<head>'
    print >> out, '<meta http-equiv="Pragma" content="no-cache">'
    print >> out, '<meta http-equiv="Expires" content="-1">'
    print >> out, '<script type="text/javascript" src="%s">' % dygraphs_url
    print >> out, '</script>'
    print >> out, '</head>'
    print >> out, '<body>'
    print >> out, '<div id="graphdiv"></div>'
    print >> out, '<script type="text/javascript">'
    print >> out, script
    print >> out, '</script>'
    print >> out, '</body>'
    print >> out, '</html>'
    return out.getvalue()

