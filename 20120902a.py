"""
Try to reproduce fig (4) from a manuscript by Nasrallah.
"""

from StringIO import StringIO
import math

import numpy as np
from scipy import linalg

import Form
import FormOut
import RUtil
from RUtil import mk_call_str
import MatrixUtil
import wfengine
import wfcompens
import multinomstate

def get_form():
    """
    @return: the body of a form
    """
    return [
            Form.CheckGroup('axis_options', 'axis options', [
                Form.CheckItem('ylogscale', 'y axis log scale'),
                ]),
            Form.Float('Nr', 'recombination parameter Nr',
                5.0, low_inclusive=0.0, high_inclusive=5.0),
            Form.ImageFormat(),
            ]

def get_form_out():
    return FormOut.Image('plot')

def rev(x):
    return tuple(reversed(list(x)))

def get_plot_array(N_diploid, Nr, theta_values, Ns_values):
    """
    Compute expected hitting times.
    Theta is 4*N*mu, and the units of time are 4*N*mu generations.
    @param N_diploid: diploid population size
    @param Nr: recombination rate
    @param theta_values: mutation rates
    @param Ns_values: selection values
    @return: arr[i][j] gives time for Ns_values[i] and theta_values[j]
    """
    # set up the state space
    k = 4
    M = multinomstate.get_sorted_states(2*N_diploid, k)
    T = multinomstate.get_inverse_map(M)
    nstates = M.shape[0]
    lmcs = wfengine.get_lmcs(M)
    # precompute rate matrices
    R_rate = wfcompens.create_recomb(M, T)
    M_rate = wfcompens.create_mutation(M, T)
    #
    r = Nr / float(N_diploid)
    #
    arr = []
    for theta in theta_values:
        mu = theta / float(4*N_diploid)
        row = []
        for Ns in Ns_values:
            s = Ns / float(N_diploid)
            lps = wfcompens.create_selection(s, M)
            S_prob = np.exp(wfengine.create_genic(lmcs, lps, M))
            M_prob = linalg.expm(theta * M_rate / float(2*2*N_diploid))
            #M_prob = linalg.expm(theta * M_rate / float(2))
            R_prob = linalg.expm(Nr * R_rate / float((2*N_diploid)**2))
            MR_prob = np.dot(M_prob, R_prob)
            P = np.dot(MR_prob, S_prob)
            # compute the stationary distribution
            v = MatrixUtil.get_stationary_distribution(P)
            # compute the transition matrix limit at time infinity
            P_inf = np.outer(np.ones_like(v), v)
            # compute the fundamental matrix Z
            Z = linalg.inv(np.eye(nstates) - (P - P_inf)) - P_inf
            # compute the hitting time from state AB to state ab.
            i = 0
            j = 3
            hitting_time_generations = (Z[j,j] - Z[i,j]) / v[j]
            hitting_time = hitting_time_generations / theta
            row.append(hitting_time)
        arr.append(row)
    return arr

def get_response_content(fs):
    # define some fixed values
    N_diploid = 4
    Nr = fs.Nr
    #plot_density = 10
    plot_density = 2
    # define some mutation rates
    theta_values = [0.001, 0.01, 0.1, 1.0]
    # define some selection coefficients to plot
    Ns_low = 0.0
    Ns_high = 3.0
    Ns_values = np.linspace(Ns_low, Ns_high, 3*plot_density + 1)
    # get the values for each h
    arr = get_plot_array(N_diploid, Nr, theta_values, Ns_values)
    # define x and y plot limits
    xlim = (Ns_low, Ns_high)
    ylim = (np.min(arr), np.max(arr))
    if fs.ylogscale:
        ylogstr = '"y"'
    else:
        ylogstr = '""'
    # colors
    colors = ['darkviolet', 'blue', 'green', 'yellow']
    # define the r script
    out = StringIO()
    print >> out, 'Ns.values <- c', str(tuple(Ns_values))
    print >> out, 'ha <- c', str(tuple(arr[0]))
    print >> out, 'hb <- c', str(tuple(arr[1]))
    print >> out, 'hc <- c', str(tuple(arr[2]))
    print >> out, 'hd <- c', str(tuple(arr[3]))
    print >> out, mk_call_str('plot', 'Ns.values', 'ha',
            type='"l"',
            xlab='"Ns"',
            ylab='"generations / theta"',
            log=ylogstr,
            main='"mean hitting time, 2N=%s"' % (2 * N_diploid),
            xlim='c' + str(xlim),
            ylim='c' + str(ylim),
            col='"%s"' % colors[0],
            )
    print >> out, mk_call_str(
            'lines', 'Ns.values', 'hb', col='"%s"' % colors[1])
    print >> out, mk_call_str(
            'lines', 'Ns.values', 'hc', col='"%s"' % colors[2])
    print >> out, mk_call_str(
            'lines', 'Ns.values', 'hd', col='"%s"' % colors[3])
    print >> out, mk_call_str(
            'legend',
            '"topleft"',
            'c' + str(tuple('theta = %s' % x for x in theta_values)),
            lty='c' + str(tuple([1]*4)),
            lwd='c' + str(tuple([2.5]*4)),
            col='c' + str(tuple(colors)),
            )
    script = out.getvalue().rstrip()
    # create the R plot image
    device_name = Form.g_imageformat_to_r_function[fs.imageformat]
    retcode, r_out, r_err, image_data = RUtil.run_plotter_no_table(
            script, device_name)
    if retcode:
        raise RUtil.RError(r_err)
    return image_data

