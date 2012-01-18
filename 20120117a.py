"""
Plot mutual information decay for a reversible continuous time Markov process.

When a reversible finite-state continuous-time
Markov process is at equilibrium,
two observations separated by time t
will be non-independent in the sense that
in a way that can be quantified by their mutual information.
This mutual information is also mathematically equivalent
to the expected log-likelihood ratio
of T=t vs. T=infinity given that the true time is T=t.
In this script we empirically investigate the asymptotic
properties of the decay of the mutual information.
The hypothesis is that the decay is
bounded above and below by an exponentially decreasing function
whose decay rate depends only on the smallest magnitude nonzero
eigenvalue of the rate matrix that defines the Markov process.
"""

from StringIO import StringIO
import math

import numpy as np
import scipy
from scipy import linalg

import Form
import FormOut
from MatrixUtil import ndot
import mrate
import divtime

def get_form():
    """
    @return: the body of a form
    """
    # define the form objects
    form_objects = [
            Form.Integer('nstates', 'number of states', 4, low=2, high=9),
            Form.Float('divtime', 'arbitrary large-ish divergence time',
                '3', low_exclusive=0)]
            #Form.Float('delta', 'vanishing time delta',
                #'0.0001', high_exclusive=0.25)]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):
    out = StringIO()
    np.set_printoptions(linewidth=200)
    # get the user defined variables
    n = fs.nstates
    t = fs.divtime
    #h = fs.delta
    # sample a random rate matrix
    v = divtime.sample_distribution(n)
    S = divtime.sample_symmetric_rate_matrix(n)
    R = divtime.to_gtr_halpern_bruno(S, v)
    # get some properties of the rate matrix
    distn = mrate.R_to_distn(R)
    spectrum = np.linalg.eigvalsh(mrate.symmetrized(R))
    #spectrum, U = np.linalg.eigh(mrate.symmetrized(R))
    #spectrum = np.linalg.eigvals(R)
    # report some information about the mutual information curve
    #mi_a = divtime.get_mutual_information(R, t - h)
    mi_b = divtime.get_mutual_information(R, t)
    #mi_c = divtime.get_mutual_information(R, t + h)
    if mi_b < 0:
        msg = 'non-positive numerically computed mutual information: %f' % mi_b
        raise ValueError(msg)
    mi_diff_b = divtime.get_mutual_information_diff(R, t)
    print >> out, 'arbitrary large-ish divergence time:'
    print >> out, t
    print >> out
    print >> out, 'randomly sampled reversible rate matrix:'
    print >> out, R
    print >> out
    print >> out, 'stationary distribution:'
    print >> out, distn
    print >> out
    print >> out, 'spectrum of the rate matrix:'
    print >> out, spectrum
    print >> out
    print >> out, 'mutual information at t = %f:' % t
    print >> out, mi_b
    print >> out
    print >> out, 'large t approximation of MI at t = %f:' % t
    print >> out, divtime.get_mutual_information_approx(R, t)
    print >> out
    print >> out, 'large t approximation of MI at t = %f (ver. 2):' % t
    print >> out, divtime.get_mutual_information_approx_b(R, t)
    print >> out
    print >> out, 'large t approximation of MI at t = %f (ver. 3):' % t
    print >> out, divtime.cute_MI_alternate(R, t)
    print >> out
    print >> out, 'mutual information diff at t = %f:' % t
    print >> out, mi_diff_b
    print >> out
    print >> out, 'large t approximation of MI diff at t = %f:' % t
    print >> out, divtime.get_mutual_information_diff_approx(R, t)
    print >> out
    print >> out, 'large t approximation of MI diff at t = %f: (ver. 2)' % t
    print >> out, divtime.get_mutual_information_diff_approx_b(R, t)
    print >> out
    print >> out, 'log of mutual information at t = %f:' % t
    print >> out, math.log(mi_b)
    print >> out
    #print >> out, 'estimated derivative',
    #print >> out, 'of log of mutual information at t = %f:' % t
    #print >> out, (math.log(mi_c) - math.log(mi_a)) / (2*h)
    #print >> out
    print >> out, 'estimated derivative of log of MI',
    print >> out, 'at t = %f:' % t
    print >> out, mi_diff_b / mi_b
    print >> out
    print >> out, 'large t approximation of derivative of log of MI',
    print >> out, 'at t = %f:' % t
    print >> out, divtime.get_mutual_information_diff_approx(R,
            t) / divtime.get_mutual_information_approx(R, t)
    print >> out
    print >> out, 'large t approximation of derivative of log of MI',
    print >> out, 'at t = %f (ver. 2):' % t
    print >> out, divtime.get_mutual_information_diff_approx_b(R,
            t) / divtime.get_mutual_information_approx_b(R, t)
    print >> out
    print >> out, 'twice the relevant eigenvalue:'
    print >> out, 2 * spectrum[-2]
    print >> out
    print >> out
    #print >> out, 'estimated derivative',
    #print >> out, 'of mutual information at t = %f:' % t
    #print >> out, (mi_c - mi_a) / (2*h)
    #print >> out
    #print >> out, '(estimated derivative of mutual information) /',
    #print >> out, '(mutual information) at t = %f:' % t
    #print >> out, (mi_c - mi_a) / (2*h*mi_b)
    #print >> out
    return out.getvalue()

