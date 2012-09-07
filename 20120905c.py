"""
Try to reproduce fig (9) from a manuscript by Nasrallah. [UNFINISHED]
"""

#TODO this needs to be changed because the variance is wrong.
# The actual variance is of a mixture model
# whose properties we know.
# we know the mixture parameter that is the probability of type 1 transition,
# and we know the variance and expectation of the number of transitions
# for type 1 and type 2 substitutions.
# the variance of the mixture can be computed using the
# principle of total variance.
#
# well now I'm not sure what is the problem
# but this does not work right now.

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
import hittingtime

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

def get_type_2_info(P):
    """
    The expected time for a type 2 event is computed as follows.
    It is the expected number of steps from AB to ab
    conditional on not entering the states AB, Ab, or aB.
    @param P: a huge transition matrix which is not modified
    @return: expectation and variance of compensatory substitution time
    """
    MatrixUtil.assert_transition_matrix(P)
    nstates = len(P)
    # define index sequences
    plain = range(4, nstates)
    forbidden = [0, 1, 2]
    target = [3]
    #
    H = hittingtime.get_conditional_transition_matrix(
            P, plain, forbidden, target)
    t = hittingtime.get_absorption_time(
            H, plain+forbidden, target)
    v = hittingtime.get_absorption_variance(
            H, plain+forbidden, target)
    return t[0], v[0]

def get_type_1_info(P):
    """
    The expected time for a type 1 event is computed as follows.
    It is the sum of two expected times.
    The first time is the expected number of steps from AB to either Ab or aB
    conditional on not entering the states AB or ab.
    The second time is the expected number of steps from Ab to ab
    conditional on not entering AB.
    Note that this formulation depends on the assumption
    that the process associated with the first
    step is equally likely to end up in Ab as in aB.
    @param P: a huge transition matrix which is not modified
    @return: expectation and variance of compensatory substitution time
    """
    MatrixUtil.assert_transition_matrix(P)
    nstates = len(P)
    # get the expected time for the first stage
    plain = range(4, nstates)
    forbidden = [0, 3]
    target = [1, 2]
    H = hittingtime.get_conditional_transition_matrix(
            P, plain, forbidden, target)
    t = hittingtime.get_absorption_time(
            H, plain+forbidden, target)
    v = hittingtime.get_absorption_variance(
            H, plain+forbidden, target)
    t_first = t[0]
    v_first = v[0]
    # get the expected time for the second stage
    plain = [1, 2] + range(4, nstates)
    forbidden = [0]
    target = [3]
    H = hittingtime.get_conditional_transition_matrix(
            P, plain, forbidden, target)
    t = hittingtime.get_absorption_time(
            H, plain+forbidden, target)
    v = hittingtime.get_absorption_variance(
            H, plain+forbidden, target)
    t_second = t[1]
    v_second = v[1]
    # return the moments of the distribution accounting for both stages
    return t_first + t_second, v_first + v_second

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
            #
            t1, v1 = get_type_1_info(P)
            t2, v2 = get_type_2_info(P)
            #
            """
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
            """
            # Get the probability of type 2.
            # This uses the stochastic complement.
            # Wait this is wrong.
            # This is the probability of a direct transition.
            X = linalg.solve(np.eye(nstates - k) - P[k:, k:], P[k:, :k])
            H = P[:k, :k] + np.dot(P[:k, k:], X)
            p_direct = H[0, 3] / (1 - H[0,0])
            # The following line is Equation (1) of the Nasrallah manuscript.
            p_t2 = (2*p_direct) / (1 + p_direct)
            p_t1 = 1 - p_t2
            """
            expectation_of_variance = p_t2*v2 + p_t1*v1
            variance_of_expectation = p_t2*p_t1*(t1 - t2)*(t1 - t2)
            pooled_variance = (
                    expectation_of_variance + variance_of_expectation)
            """
            #
            # Just do a simulation,
            # assuming that the wait times are normally distributed.
            nsamples = 500
            n1 = np.random.binomial(nsamples, p_t1)
            n2 = nsamples - n1
            X1 = np.random.normal(t1, math.sqrt(v1), n1)
            X2 = np.random.normal(t2, math.sqrt(v2), n2)
            X_pooled = np.hstack((X1, X2))
            x = np.mean(X1) - np.mean(X2)
            s_pooled = math.sqrt(np.var(X_pooled) / nsamples)
            t_statistic = x / s_pooled
            row.append(t_statistic)
            #
            #x = (t1 - t2) / math.sqrt(variance / 200.0)
            #x = (t1 - t2) / math.sqrt((v1 + v2) / 200.0)
            #x = (t1 - t2) / math.sqrt(pooled_variance)
            #x = (t1 - t2)
            #row.append(math.log(t1) - math.log(t2))
            #row.append(x)
            #row.append(v2)
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
            type='"p"',
            xlab='"Ns"',
            ylab=ylab,
            log=ylogstr,
            main='"theta=%s ; 2N=%s"' % (theta, N_hap),
            xlim='c' + str(xlim),
            ylim='c' + str(ylim),
            col='"%s"' % colors[0],
            )
    print >> out, mk_call_str(
            'points', 'Ns.values', 'hb', col='"%s"' % colors[1])
    print >> out, mk_call_str(
            'legend',
            '"topright"',
            'c' + str(tuple('%s' % x for x in Nr_values)),
            title='"Nr"',
            lty='c' + str(tuple([1]*2)),
            lwd='c' + str(tuple([2.5]*2)),
            col='c' + str(tuple(colors)),
            )
    return out.getvalue().rstrip()

def get_response_content(fs):
    # define some fixed values
    N_diploid = 8
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
    Ns_high = 2.5
    Ns_values = np.linspace(Ns_low, Ns_high, 3*plot_density + 1)
    # get the values for each h
    arr = get_plot_array(N_diploid, theta, Nr_values, Ns_values)
    #ylab='"log(Type1 / Type2)"'
    ylab='"normalized time (Type1 - Type2)"'
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

