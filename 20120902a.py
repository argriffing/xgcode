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
            Form.CheckGroup('misc_options', 'misc options', [
                Form.CheckItem('ylogscale', 'y axis log scale'),
                Form.CheckItem('scale_to_2N_200', 'scale to match 2N=200'),
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
    # precompute a recombination probability matrix
    R_prob = linalg.expm(Nr * R_rate / float((2*N_diploid)**2))
    #
    arr = []
    for theta in theta_values:
        # Compute the expected number of mutation events per generation.
        mu = theta / 2
        # Precompute the mutation matrix
        # and the product of mutation and recombination.
        M_prob = linalg.expm(mu * M_rate / float(2*2*N_diploid))
        MR_prob = np.dot(M_prob, R_prob)
        #
        row = []
        for Ns in Ns_values:
            s = Ns / float(N_diploid)
            lps = wfcompens.create_selection(s, M)
            S_prob = np.exp(wfengine.create_genic(lmcs, lps, M))
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
            hitting_time = hitting_time_generations * theta
            row.append(hitting_time)
        arr.append(row)
    return arr

def get_plot(
        side, Nr, arr, theta_values, Ns_values, xlim, ylim, ylogstr, ylab):
    """
    @param side: 'left' or 'right'
    """
    colors = ['blue', 'pink', 'green', 'red']
    out = StringIO()
    print >> out, 'Ns.values <- c', str(tuple(Ns_values))
    print >> out, 'ha <- c', str(tuple(arr[0]))
    print >> out, 'hb <- c', str(tuple(arr[1]))
    print >> out, 'hc <- c', str(tuple(arr[2]))
    print >> out, 'hd <- c', str(tuple(arr[3]))
    print >> out, mk_call_str('plot', 'Ns.values', 'ha',
            type='"l"',
            xlab='"Ns"',
            ylab=ylab,
            log=ylogstr,
            main='"Nr=%s"' % Nr,
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
    if side == 'left':
        print >> out, mk_call_str(
                'legend',
                '"topleft"',
                'c' + str(tuple('%s' % x for x in theta_values)),
                title='"theta"',
                lty='c' + str(tuple([1]*4)),
                lwd='c' + str(tuple([2.5]*4)),
                col='c' + str(tuple(colors)),
                )
    return out.getvalue().rstrip()

def get_response_content(fs):
    # define some fixed values
    N_diploid = 10
    N_hap = 2 * N_diploid
    #Nr = fs.Nr
    plot_density = 2
    # define some mutation rates
    theta_values = [0.001, 0.01, 0.1, 1.0]
    # define some selection coefficients to plot
    Ns_low = 0.0
    Ns_high = 3.0
    Ns_values = np.linspace(Ns_low, Ns_high, 3*plot_density + 1)
    # get the values for each h
    Nr_values = (0, 5)
    arr_0 = get_plot_array(
            N_diploid, Nr_values[0], theta_values, Ns_values)
    arr_1 = get_plot_array(
            N_diploid, Nr_values[1], theta_values, Ns_values)
    if fs.scale_to_2N_200:
        arr_0 = (200 / float(N_hap)) * np.array(arr_0)
        arr_1 = (200 / float(N_hap)) * np.array(arr_1)
        ylab = '"generations * theta * (200 / 2N)"'
    else:
        ylab='"generations * theta"'
    # define x and y plot limits
    xlim = (Ns_low, Ns_high)
    ylim = (np.min((arr_0, arr_1)), np.max((arr_0, arr_1)))
    if fs.ylogscale:
        ylogstr = '"y"'
    else:
        ylogstr = '""'
    # http://sphaerula.com/legacy/R/multiplePlotFigure.html
    out = StringIO()
    print >> out, mk_call_str(
            'par',
            mfrow='c(1,2)',
            oma='c(0,0,2,0)',
            )
    print >> out, get_plot(
            'left', Nr_values[0], arr_0, theta_values, Ns_values,
            xlim, ylim, ylogstr, ylab)
    print >> out, get_plot(
            'right', Nr_values[1], arr_1, theta_values, Ns_values,
            xlim, ylim, ylogstr, '""')
    print >> out, mk_call_str(
            'title',
            '"mean hitting time, 2N=%s"' % N_hap,
            outer='TRUE',
            )
    script = out.getvalue().rstrip()
    # create the R plot image
    device_name = Form.g_imageformat_to_r_function[fs.imageformat]
    retcode, r_out, r_err, image_data = RUtil.run_plotter_no_table(
            script, device_name)
    if retcode:
        raise RUtil.RError(r_err)
    return image_data


def main():
    """
    Check some stuff that takes longer to compute.
    """
    pass

if __name__ == '__main__':
    main()

