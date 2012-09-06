"""
Try to reproduce fig (9) from a manuscript by Nasrallah.
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
            Form.RadioGroup('theta_choice', 'choice of theta', [
                Form.RadioItem('theta_1em0', '1', True),
                Form.RadioItem('theta_1em1', '0.1'),
                Form.RadioItem('theta_1em2', '0.01'),
                ]),
            Form.ImageFormat(),
            ]

def get_form_out():
    return FormOut.Image('plot')

def get_absorption_time(P_absorbing):
    """
    This function assumes that exactly one state is absorbing.
    It returns the time to absorption from state 0 to state -1.
    The final state is assumed to be the absorbing state
    @param P_absorbing: an absorbing Markov chain
    @return: time to absorption from state 0
    """
    nstates = P_absorbing.shape[0]
    Q = P_absorbing[:-1, :-1]
    t = linalg.solve(np.eye(nstates-1) - Q, np.ones(nstates-1))
    return t[0]

def get_type_2_absorption_time(P):
    """
    The indices of the transition matrix are ordered as follows.
    The last three states are fixed Ab, fixed aB, and fixed ab.
    In type 2 events we have the following conditions.
    Nothing transitions to fixed AB except for fixed AB itself.
    Nothing transitions to fixed Ab or fixed aB.
    Fixed ab is an absorbing state.
    @param P: a mutable copy of the huge transition matrix
    @return: expected time to absorption from fixed state AB
    """
    nstates = P.shape[0]
    # do not allow transitions to the AB state except from itself
    P[1:, 0] = 0
    # do not allow transitions to the low fitness fixed states
    P[:, -3:-1] = 0
    # force the last three states to be absorbing
    P[-3:] = 0
    P[-3:, -3:] = np.eye(3)
    # normalize the rows
    v = P.sum(axis=1)
    P /= v[:, np.newaxis]
    # Compute the time to absorption using standard notation.
    Q = P[:-3, :-3]
    c = np.ones(nstates-3)
    t = linalg.solve(np.eye(nstates-3) - Q, c)
    return t[0]

def get_type_1_absorption_time(P):
    """
    Ordering of indices is the same as in type 2.
    In type 1 events we have the following conditions.
    Nothing transitions to fixed AB except for fixed AB itself.
    Fixed ab is an absorbing state.
    @param P: a mutable copy of the huge transition matrix
    @return: expected time to absorption from fixed state AB
    """
    nstates = P.shape[0]
    # do not allow transitions to the AB state except from itself
    P[1:, 0] = 0
    # force the ab state to be absorbing
    P[-1] = 0
    P[-1, -1] = 1
    # normalize the rows
    v = P.sum(axis=1)
    P /= v[:, np.newaxis]
    # Compute the time to absorption using standard notation.
    Q = P[:-1, :-1]
    c = np.ones(nstates-1)
    t = linalg.solve(np.eye(nstates-1) - Q, c)
    return t[0]

def get_plot_array(N_diploid, theta, Nr_values, Ns_values):
    """
    @param N_diploid: diploid population size
    @param theta: mutation rate
    @param Nr_values: recombination rates
    @param Ns_values: selection values
    @return: arr[i][j] gives time for Ns_values[i] and theta_values[j]
    """
    # define the haplotypes
    AB, Ab, aB, ab = 0, 1, 2, 3
    # initialize the state space
    N_hap = 2 * N_diploid
    k = 4
    M = multinomstate.get_sorted_states(N_hap, k)
    nstates = M.shape[0]
    # Swap the fixed states to the end
    # to more closely match standard forms.
    M[:4], M[-4:] = M[-4:], M[:4]
    # compute the inverse map
    T = multinomstate.get_inverse_map(M)
    #
    lmcs = wfengine.get_lmcs(M)
    # precompute rate matrices
    R_rate = wfcompens.create_recomb(M, T)
    M_rate = wfcompens.create_mutation(M, T)
    # Compute the expected number of mutation events per generation.
    mu = theta / 2
    # Precompute the mutation matrix
    M_prob = linalg.expm(mu * M_rate / float(2*2*N_diploid))
    #
    arr = []
    for Nr in Nr_values:
        # precompute a recombination probability matrix
        R_prob = linalg.expm(Nr * R_rate / float((2*N_diploid)**2))
        # precompute the product of mutation and recombination.
        MR_prob = np.dot(M_prob, R_prob)
        #
        row = []
        for Ns in Ns_values:
            s = Ns / float(N_diploid)
            lps = wfcompens.create_selection(s, M)
            S_prob = np.exp(wfengine.create_genic(lmcs, lps, M))
            P = np.dot(MR_prob, S_prob)
            # What is the distribution over next fixed states
            # from the current state?
            # This question can be answered
            # by hacking with transience and absorption.
            Q = P[:-k, :-k]
            R = P[:-k, -k:]
            B = linalg.solve(np.eye(nstates-k) - Q, R)
            # At this point B is the matrix whose nstates-k rows give
            # distributions over the k fixed states.
            # Next construct the transition matrix that is conditional
            # upon first hitting the ab fixed state.
            w = np.zeros(nstates)
            w[:-k] = R[:, -1]
            w[-k:] = np.array([0, 0, 0, 1])
            P_t2 = P * w
            # normalize after scaling by the weights
            v = P_t2.sum(axis=1)
            P_t2 /= v[:, np.newaxis]
            # Get the hitting time from state AB to state ab.
            # Because of the conditioning, this should be the same
            # as the expected time to reach state ab given that state ab
            # is the first reached fixed state.
            # Note that this means that the first step is away from AB.
            # Or actually we can just use expected time to absorption.
            Q = P_t2[:-1, :-1]
            c = np.ones(nstates-1)
            t = linalg.lstsq(np.eye(nstates-1) - Q, c)
            t2 = t[-4]
            # Now do type 1 events.
            w = np.zeros(nstates)
            w[:-k] = 1 - R[:, 0]
            w[-k:] = np.array([0, 0, 0, 1])
            P_t2 = P * w
            #
            row.append(math.log(t1) - math.log(t2))
        arr.append(row)
    return arr

def get_plot(
        N_hap, theta, arr, Nr_values, Ns_values, xlim, ylim, ylogstr, ylab):
    colors = ['blue', 'red']
    out = StringIO()
    print >> out, 'Ns.values <- c', str(tuple(Ns_values))
    print >> out, 'ha <- c', str(tuple(arr[0]))
    print >> out, 'hb <- c', str(tuple(arr[1]))
    print >> out, mk_call_str('plot', 'Ns.values', 'ha',
            type='"l"',
            xlab='"Ns"',
            ylab=ylab,
            log=ylogstr,
            main='"theta=%s ; 2N=%s"' % (theta, N_hap),
            xlim='c' + str(xlim),
            ylim='c' + str(ylim),
            col='"%s"' % colors[0],
            )
    print >> out, mk_call_str(
            'lines', 'Ns.values', 'hb', col='"%s"' % colors[1])
    print >> out, mk_call_str(
            'legend',
            '"topleft"',
            'c' + str(tuple('%s' % x for x in Nr_values)),
            title='"Nr"',
            lty='c' + str(tuple([1]*2)),
            lwd='c' + str(tuple([2.5]*2)),
            col='c' + str(tuple(colors)),
            )
    return out.getvalue().rstrip()

def get_response_content(fs):
    # define some fixed values
    N_diploid = 5
    N_hap = 2 * N_diploid
    plot_density = 8
    # get the user-defined theta
    if fs.theta_1em0:
        theta = 1.0
    elif fs.theta_1em1:
        theta = 0.1
    elif fs.theta_1em2:
        theta = 0.01
    # define some mutation rates
    Nr_values = [0.0, 5.0]
    # define some selection coefficients to plot
    Ns_low = 0.0
    Ns_high = 3.0
    Ns_values = np.linspace(Ns_low, Ns_high, 3*plot_density + 1)
    # get the values for each h
    arr = get_plot_array(N_diploid, theta, Nr_values, Ns_values)
    ylab='"log(Type1 / Type2)"'
    # define x and y plot limits
    xlim = (Ns_low, Ns_high)
    ylim = (np.min(arr), np.max(arr))
    ylogstr = '""'
    out = StringIO()
    print >> out, get_plot(
            N_hap, theta, arr, Nr_values, Ns_values,
            xlim, ylim, ylogstr, ylab)
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

