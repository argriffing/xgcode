r"""
Compare time information with and without selection.

The way that the GTR mutation matrix is broken into the
strict lower diagonal exchangeability component and the
stationary distribution vector
(here unnormalized for the convenience of not having to manually
re-normalize the distribution every time you want to change it)
is taken from the wikipedia article on models of DNA evolution.
In particular, \( \mu_{xy} = S_{xy} \pi_y \).
The way that the mutation-selection balance stationary distribution
is used to construct the mutation-selection balance rate matrix
given the GTR mutation matrix is taken from
the paper by Sang-Chul Choi et al.
"""

from StringIO import StringIO
import argparse

import numpy as np
import scipy
from scipy import linalg

import Form
import FormOut
import iterutils
import ctmcmi
import mrate
import divtime

g_default_exch = """
1
1 1
1 1 1
"""

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('lowtri',
                'strictly lower triangular mut exch',
                g_default_exch.strip()),
            Form.Sequence('mutweights',
                'unnormalized mut stationary distn',
                ('1', '1', '1', '1')),
            Form.Sequence('mutselweights',
                'unnormalized mut-sel stationary distn',
                ('1', '1', '0.01', '0.01')),
            Form.Float('t', 'time', '0.1', low_exclusive=0),
            Form.RadioGroup('infotype', 'information variant', [
                Form.RadioItem('info_fis', 'Fisher information', True),
                Form.RadioItem('info_mut', 'mutual information')])]
    return form_objects

def get_form_out():
    return FormOut.Report('summary')

def get_input_matrices(fs):
    """
    @return: M, R
    """
    # get the positive strict lower triangular part of the S matrix
    L = []
    for i, line in enumerate(iterutils.stripped_lines(fs.lowtri.splitlines())):
        values = line.split()
        if len(values) != i + 1:
            msg = 'expected %d values on line "%s"' % (
                    i+1, line)
            raise ValueError(msg)
        vs = [float(v) for v in values]
        if any(x<0 for x in vs):
            raise ValueError('exchangeabilities must be nonnegative')
        L.append(vs)
    # get the mut and mutsel weights
    mut_weights = [float(v) for v in fs.mutweights]
    mutsel_weights = [float(v) for v in fs.mutselweights]
    if any(x<=0 for x in mut_weights + mutsel_weights):
        raise ValueError('stationary weights must be positive')
    # normalize weights to distributions
    mut_distn = [v / sum(mut_weights) for v in mut_weights]
    mutsel_distn = [v / sum(mutsel_weights) for v in mutsel_weights]
    # get the time
    t = fs.t
    # get the exchangeability matrix
    nstates = len(L) + 1
    S = np.zeros((nstates, nstates))
    for i, row in enumerate(L):
        for j, v in enumerate(row):
            S[i+1, j] = v
            S[j, i+1] = v
    # check the state space sizes implied by the inputs
    if len(set(len(x) for x in (S, mut_weights, mutsel_weights))) != 1:
        msg = 'the inputs do not agree on the state space size'
        raise ValueError(msg)
    # check for sufficient number of states
    if nstates < 2:
        msg = 'at least two states are required'
        raise ValueError(msg)
    # get the mutation rate matrix
    M = np.zeros((nstates, nstates))
    for i in range(nstates):
        for j in range(nstates):
            M[i, j] = S[i, j] * mut_distn[j]
    M -= np.diag(np.sum(M, axis=1))
    # get the mutation selection balance rate matrix
    R = ctmcmi.to_gtr_halpern_bruno(M, mutsel_distn)
    # check the stationary distributions
    mut_distn_observed = mrate.R_to_distn(M)
    if not np.allclose(mut_distn_observed, mut_distn):
        msg = 'internal mut stationary distribution computation error'
        raise ValueError(msg)
    mutsel_distn_observed = mrate.R_to_distn(R)
    if not np.allclose(mutsel_distn_observed, mutsel_distn):
        msg = 'internal mut-sel stationary distribution computation error'
        raise ValueError(msg)
    # return the values
    return M, R

def get_response_content(fs):
    M, R = get_input_matrices(fs)
    M_v = mrate.R_to_distn(M)
    R_v = mrate.R_to_distn(R)
    t = fs.t
    if fs.info_mut:
        mut_information = ctmcmi.get_mutual_information(M, t)
        mutsel_information = ctmcmi.get_mutual_information(R, t)
    elif fs.info_fis:
        mut_information = divtime.get_fisher_information(M, t)
        mutsel_information = divtime.get_fisher_information(R, t)
    # print a summary
    out = StringIO()
    print >> out, 'mutation rate matrix:'
    print >> out, M
    print >> out
    print >> out, 'mutation process stationary distribution:'
    print >> out, M_v
    print >> out
    print >> out, 'mutation-selection balance rate matrix:'
    print >> out, R
    print >> out
    print >> out, 'mutation-selection balance stationary distribution:'
    print >> out, R_v
    print >> out
    print >> out, 'mutation process expected rate:'
    print >> out, mrate.Q_to_expected_rate(M)
    print >> out
    print >> out, 'mutation-selection balance expected rate:'
    print >> out, mrate.Q_to_expected_rate(R)
    print >> out
    print >> out, 'mutation process information (t = %s):' % t
    print >> out, mut_information
    print >> out
    print >> out, 'mutation-selection balance information (t = %s):' % t
    print >> out, mutsel_information
    print >> out
    print >> out, 'zone:'
    if mut_information > mutsel_information:
        print >> out, 'pure mutation gives more information'
    else:
        print >> out, 'mutation-selection balance gives more information'
    return out.getvalue()

