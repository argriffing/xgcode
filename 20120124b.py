"""
Plot mutual information using the dygraph javascript.
"""


from StringIO import StringIO
import math
import random

import numpy as np
import scipy
from scipy import stats

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
    return FormOut.Html('timelines')

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
    # add the mutual information sign proportion
    statistics.append(sum(1 for x in mi_signs if x == 1) / float(nsels))
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

def get_response_content(fs):
    # define remote mathjax code
    mathjax_js = 'cdn.mathjax.org/mathjax/latest/MathJax.js'
    mathjax_params = 'config=TeX-AMS-MML_HTMLorMML'
    mathjax_url = 'http://' + mathjax_js + '?' + mathjax_params
    # define remote dygraph code
    dygraphs_url = 'http://dygraphs.com/dygraph-combined.js'
    time_stats = get_time_stats(fs)
    data = zip(*time_stats)
    # create an html file that uses the javascript link
    out = StringIO()
    print >> out, '<html>'
    print >> out, '<head>'
    print >> out, '<meta http-equiv="Pragma" content="no-cache">'
    print >> out, '<meta http-equiv="Expires" content="-1">'
    print >> out, '<script type="text/javascript" src="%s">' % mathjax_url
    print >> out, '</script>'
    print >> out, '<script type="text/javascript" src="%s">' % dygraphs_url
    print >> out, '</script>'
    print >> out, '</head>'
    print >> out, '<body>'
    #
    print >> out, '<p>'
    print >> out, 'The following plot shows the mutual information'
    print >> out, 'between observations'
    print >> out, 'separated by increasing amounts of time.'
    print >> out, 'The observations are from reversible finite-state'
    print >> out, 'continuous-time Markov processes.'
    print >> out, 'The mutual information has been computed for'
    print >> out, 'a mutation process'
    print >> out, 'and for a sample of %d corresponding' % fs.nselections
    print >> out, 'mutation-selection balance processes.'
    print >> out, 'For the mutation-selection balance processes,'
    print >> out, 'only the mean and the range of the mutual information'
    print >> out, 'computations are shown for each'
    print >> out, 'separation time.'
    print >> out, get_dygraph_html(get_script_a, 'gdiv1', data)
    print >> out, '</p>'
    #
    print >> out, '<p>'
    print >> out, 'The following plot shows the proportion'
    print >> out, 'of mutation-selection processes'
    print >> out, 'for which the mutual information is greater'
    print >> out, 'than that of the pure mutation process.'
    print >> out, get_dygraph_html(get_script_b, 'gdiv2', data)
    print >> out, '</p>'
    #
    print >> out, '<p>'
    print >> out, 'The next three plots compare'
    print >> out, 'five different functions to the mutual information.'
    print >> out, r'Let \( \pi_i \) be the stationary probability'
    print >> out, r'of state \( i \) under the mutation-selection process.'
    print >> out, r'Let \( P(t)_{i,j} \) be the conditional probability'
    print >> out, r'of ending up in state \( j \) after time \( t \)'
    print >> out, r'given that the starting state is \( i \),'
    print >> out, r'under the mutation-selection balance process.'
    print >> out, r'Let \( \lambda_2 \) be the "second eigenvalue"'
    print >> out, 'of the mutation-selection rate matrix;'
    print >> out, 'if there is a distinct eigenvalue equal to zero'
    print >> out, 'and all of the other eigenvalues are negative,'
    print >> out, 'then the "second eigenvalue" is the negative eigenvalue'
    print >> out, 'closest to zero.'
    print >> out, r'Let \( r \) be the expected rate'
    print >> out, 'of the mutation-selection process.'
    print >> out, 'The five functions are as follows:'
    print >> out, '<ul>'
    print >> out, r'<li> \( E_i \left['
    print >> out, r'     P(t)_{i,i}'
    print >> out, r'     \log \left( \frac{1}{\pi_i} \right)'
    print >> out, r'     \right] \) </li>'
    print >> out, r'<li> \( E_i \left['
    print >> out, r'     P(t)_{i,i}'
    print >> out, r'     \log \left( \frac{P(t)_{i,i}}{\pi_i} \right)'
    print >> out, r'     \right] \) </li>'
    print >> out, r'<li> \( E_{i,j} \left['
    print >> out, r'     \frac{1}{2}'
    print >> out, r'     \left( \frac{P(t)_{i,j}}{\pi_j} - 1 \right)'
    print >> out, r'     \right] \) </li>'
    print >> out, r'<li> \( e^{2 \lambda_2 t} \) </li>'
    print >> out, r'<li> \( e^{- r t} \) </li>'
    print >> out, '</ul>'
    print >> out, 'The idea is these functions are in some ways simpler'
    print >> out, 'than the mutual information,'
    print >> out, 'which is itself written as follows:'
    print >> out, r'\[ E_{i,j} \left['
    print >> out, r'   \log \left('
    print >> out, r'   \frac{P(t)_{i,j}}{\pi_j} \right)'
    print >> out, r'   \right] \]'
    print >> out, 'If these simpler functions share properties'
    print >> out, 'with the mutual information,'
    print >> out, 'then these properties of mutual information'
    print >> out, 'can be studied through these simpler functions.'
    print >> out, '<p>'
    #
    print >> out, '</p>'
    print >> out, 'The following plot shows the correlations'
    print >> out, 'of the functions with the mutual information.'
    print >> out, get_dygraph_html(get_script_c, 'gdiv3', data)
    print >> out, '</p>'
    #
    print >> out, '<p>'
    print >> out, 'The following plot shows the proportion of samples'
    print >> out, 'for which the given function differs from the'
    print >> out, 'mutual information of the mutation process'
    print >> out, 'in the same direction (positive vs. negative)'
    print >> out, 'as does the mutual information of the'
    print >> out, 'mutuation-selection process.'
    print >> out, get_dygraph_html(get_script_d, 'gdiv4', data)
    print >> out, '</p>'
    #
    print >> out, '<p>'
    print >> out, 'The following plot shows the "informativeness"'
    print >> out, 'of the given function with respect to the'
    print >> out, 'direction of mutual information change.'
    print >> out, get_dygraph_html(get_script_e, 'gdiv5', data)
    print >> out, '</p>'
    #
    print >> out, '</body>'
    print >> out, '</html>'
    return out.getvalue()

