"""
Check Wright-Fisher compensatory time expectations and variances.

Compare estimates obtained by forward simulation
to estimates obtained through numerical solutions of systems of equations
that approximate the forward simulation dynamics.
Two different compensatory substitution pathways are considered.
The model is of two alleles at two loci, so a four-haplotype model.
Mutation, recombination, and selection+drift occur between each generation
in that order.
Selection acts directly on the haploid population,
in a way that corresponds to multiplicative selection in diploids.
The two different compensatory substitution pathways
from the high-fitness state AB to the high-fitness state ab are distinguished
by whether the most recently fixed haplotype before first reaching state ab
is AB or is a low-fitness haplotype.
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
import wfbckcompens

def get_form():
    """
    @return: the body of a form
    """
    return [
            Form.Integer('N_diploid', 'diploid population size N', 6,
                low=1, high=10),
            Form.Float('theta', 'mutation rate 4*N*mu', '1.0',
                low_inclusive=0),
            Form.Float('Nr', 'recombination rate N*r', '5.0',
                low_inclusive=0),
            Form.Float('Ns', 'selection N*s', '1.0',
                low_inclusive=0),
            ]

def get_form_out():
    return FormOut.Report()

def get_type_1_info():
    """
    Get the expectation and variance for each of the types.
    """
    pass

def get_plot_array(N_diploid, Nr, theta_values, Ns_values):
    """
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
            # compute the stochastic complement
            X = linalg.solve(np.eye(nstates - k) - P[k:, k:], P[k:, :k])
            H = P[:k, :k] + np.dot(P[:k, k:], X)
            # condition on not looping
            np.fill_diagonal(H, 0)
            v = np.sum(H, axis=1)
            H /= v[:, np.newaxis]
            # let ab be an absorbing state
            # and compute the expected number of returns to AB
            Q = H[:3, :3]
            I = np.eye(3)
            N = linalg.inv(I - Q)
            #
            row.append(N[0, 0] - 1)
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
    N_diploid = 6
    N_hap = 2 * N_diploid
    plot_density = 8
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
    ylab='"expected returns to AB"'
    # define x and y plot limits
    xlim = (Ns_low, Ns_high)
    ylim = (np.min((arr_0, arr_1)), np.max((arr_0, arr_1)))
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
            '"expected number of returns to AB, 2N=%s"' % N_hap,
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

