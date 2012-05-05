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
import math

import numpy as np
import scipy
from scipy import linalg

import Form
import FormOut
import ctmcmi
import mrate
import divtime
import MatrixUtil

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.Sequence('lowtri',
                'strictly lower triangular mut exch',
                ('1', '1 1', '1 1 1')),
            Form.Sequence('mutweights',
                'unnormalized mut stationary distn',
                ('1', '1', '1', '1')),
            Form.Float('mutscale',
                'extra mutation process scaling factor',
                '1', low_exclusive=0),
            Form.Sequence('mutselweights',
                'unnormalized mut-sel stationary distn',
                ('1', '1', '0.01', '0.01')),
            Form.Float('t', 'time', '0.1', low_exclusive=0),
            Form.RadioGroup('infotype', 'information variant', [
                Form.RadioItem('info_fis', 'Fisher information', True),
                Form.RadioItem('info_mut', 'mutual information')])]
    return form_objects

def get_presets():
    presets = [
            Form.Preset(
                'a 4-state-to-2-state example with t = 0.5',
                {
                    'lowtri' : ('1', '1 1', '1 1 1'),
                    'mutweights' : ('1', '1', '1', '1'),
                    'mutscale' : '1.333333333333',
                    'mutselweights' : ('1', '1', '0.0001', '0.0001'),
                    't' : '0.5',
                    'infotype' : 'info_mut'}),
            Form.Preset(
                'a 4-state-to-2-state example with t = 1.5',
                {
                    'lowtri' : ('1', '1 1', '1 1 1'),
                    'mutweights' : ('1', '1', '1', '1'),
                    'mutscale' : '1.333333333333',
                    'mutselweights' : ('1', '1', '0.0001', '0.0001'),
                    't' : '1.5',
                    'infotype' : 'info_mut'}),
            Form.Preset(
                'four state bottleneck',
                {
                    'lowtri' : ('1', '0 1', '1 0 1'),
                    'mutweights' : ('1', '1', '1', '1'),
                    'mutscale' : '1',
                    'mutselweights' : ('1', '0.01', '1', '0.01'),
                    't' : '0.1',
                    'infotype' : 'info_fis'}),
            Form.Preset(
                'two state, mut has greater entropy',
                {
                    'lowtri' : ('1',),
                    'mutweights' : ('0.5', '0.5'),
                    'mutscale' : '1',
                    'mutselweights' : ('0.25', '0.75'),
                    't' : '0.1',
                    'infotype' : 'info_fis'}),
            Form.Preset(
                'two state, mut-sel has greater entropy',
                {
                    'lowtri' : ('1',),
                    'mutweights' : ('0.25', '0.75'),
                    'mutscale' : '1',
                    'mutselweights' : ('0.5', '0.5'),
                    't' : '0.1',
                    'infotype' : 'info_fis'})]
    return presets

def get_form_out():
    return FormOut.Report('summary')

def _heuristic_helper(m_sign):
    """
    Return a message about which selection is predicted to be more informative.
    That is, whether pure mutation or mutation-selection balance
    is predicted by the heuristic to be more informative for divergence time.
    @param m_sign: 1 when the heuristic suggests pure mutation
    @return: a single line message
    """
    suffix = 'is predicted to be more informative'
    if m_sign == 1:
        s = '* pure mutation ' + suffix
    elif m_sign == -1:
        s = '* the balance of mutation and selection ' + suffix
    else:
        s = '  this heuristic is uninformative'
    return s

def get_heuristics(M, R):
    """
    Return a multiline string with some heuristics.
    The heuristics are independendent of time and of the information variant.
    Greater stationary distribution shannon entropy suggests less saturation.
    Greater stationary distribution logical entropy suggests less saturation.
    Greater expected rate suggests more saturation.
    Greater spectral rate suggests more saturation.
    @param M: mutation rate matrix
    @param R: mutation-selection balance rate matrix
    @return: multiline string
    """
    # get the stationary distributions
    M_v = mrate.R_to_distn(M)
    R_v = mrate.R_to_distn(R)
    # check a different way to get the stationary distribution just for fun
    M_v_nonspectral = mrate.R_to_distn_nonspectral(M)
    R_v_nonspectral = mrate.R_to_distn_nonspectral(R)
    if not np.allclose(M_v, M_v_nonspectral):
        raise ValueError('internal stationary distribution calculation error')
    if not np.allclose(R_v, R_v_nonspectral):
        raise ValueError('internal stationary distribution calculation error')
    # compute the shannon entropy of the matrices
    M_shannon_entropy = -sum(p * math.log(p) for p in M_v)
    R_shannon_entropy = -sum(p * math.log(p) for p in R_v)
    shannon_entropy_sign = np.sign(M_shannon_entropy - R_shannon_entropy)
    # compute the logical entropy of the matrices
    M_logical_entropy = 1 - sum(p * p for p in M_v)
    R_logical_entropy = 1 - sum(p * p for p in R_v)
    logical_entropy_sign = np.sign(M_logical_entropy - R_logical_entropy)
    # compute the expected rate
    M_expected_rate = mrate.Q_to_expected_rate(M)
    R_expected_rate = mrate.Q_to_expected_rate(R)
    expected_rate_sign = np.sign(R_expected_rate - M_expected_rate)
    # compute the spectral rate
    M_spectral_rate = 1 / mrate.R_to_relaxation_time(M)
    R_spectral_rate = 1 / mrate.R_to_relaxation_time(R)
    spectral_rate_sign = np.sign(R_spectral_rate - M_spectral_rate)
    # report the heuristics
    out = StringIO()
    print >> out, 'Greater Shannon entropy of the stationary distribution',
    print >> out, 'suggests more information about divergence time.'
    print >> out, _heuristic_helper(shannon_entropy_sign)
    print >> out
    print >> out, 'Greater logical entropy of the stationary distribution',
    print >> out, 'suggests more information about divergence time.'
    print >> out, _heuristic_helper(logical_entropy_sign)
    print >> out
    print >> out, 'Smaller expected rate',
    print >> out, 'suggests more information about divergence time.'
    print >> out, _heuristic_helper(expected_rate_sign)
    print >> out
    print >> out, 'Smaller spectral rate',
    print >> out, 'suggests more information about divergence time.'
    print >> out, _heuristic_helper(spectral_rate_sign)
    print >> out
    return out.getvalue().strip()

def get_response_content(fs):
    M, R = get_input_matrices(fs)
    M_v = mrate.R_to_distn(M)
    R_v = mrate.R_to_distn(R)
    t = fs.t
    mi_mut = ctmcmi.get_mutual_information(M, t)
    mi_bal = ctmcmi.get_mutual_information(R, t)
    fi_mut = divtime.get_fisher_information(M, t)
    fi_bal = divtime.get_fisher_information(R, t)
    if fs.info_mut:
        information_sign = np.sign(mi_mut - mi_bal)
    elif fs.info_fis:
        information_sign = np.sign(fi_mut - fi_bal)
    out = StringIO()
    print >> out, 'Explicitly computed answer',
    print >> out, '(not a heuristic but may be numerically imprecise):'
    if information_sign == 1:
        print >> out, '* pure mutation',
        print >> out, 'is more informative'
    elif information_sign == -1:
        print >> out, '* the balance of mutation and selection',
        print >> out, 'is more informative'
    else:
        print >> out, '  the information contents of the two processes',
        print >> out, 'are numerically indistinguishable'
    print >> out
    print >> out
    print >> out, 'Heuristics without regard to time or to the selected',
    print >> out, 'information variant (Fisher vs. mutual information):'
    print >> out
    print >> out, get_heuristics(M, R)
    print >> out
    print >> out
    print >> out, 'Input summary:'
    print >> out
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
    print >> out
    print >> out, 'The following information calculations',
    print >> out, 'depend on t = %s:' % t
    print >> out
    print >> out, 'log(ratio(E(L))) for pure mutation:'
    print >> out, ctmcmi.get_ll_ratio_wrong(M, t)
    print >> out
    print >> out, 'log(ratio(E(L))) for mut-sel balance:'
    print >> out, ctmcmi.get_ll_ratio_wrong(R, t)
    print >> out
    print >> out, 'mutual information for pure mutation:'
    print >> out, mi_mut
    print >> out
    print >> out, 'mutual information for mut-sel balance:'
    print >> out, mi_bal
    print >> out
    print >> out, 'pinsker lower bound mi for pure mutation:'
    print >> out, ctmcmi.get_pinsker_lower_bound_mi(M, t)
    print >> out
    print >> out, 'pinsker lower bound mi for mut-sel balance:'
    print >> out, ctmcmi.get_pinsker_lower_bound_mi(R, t)
    print >> out
    print >> out, 'row based pinsker lower bound mi for pure mutation:'
    print >> out, ctmcmi.get_row_based_plb_mi(M, t)
    print >> out
    print >> out, 'row based pinsker lower bound mi for mut-sel balance:'
    print >> out, ctmcmi.get_row_based_plb_mi(R, t)
    print >> out
    print >> out, 'row based hellinger lower bound mi for pure mutation:'
    print >> out, ctmcmi.get_row_based_hellinger_lb_mi(M, t)
    print >> out
    print >> out, 'row based hellinger lower bound mi for mut-sel balance:'
    print >> out, ctmcmi.get_row_based_hellinger_lb_mi(R, t)
    print >> out
    print >> out, 'Fisher information for pure mutation:'
    print >> out, fi_mut
    print >> out
    print >> out, 'Fisher information for mut-sel balance:'
    print >> out, fi_bal
    print >> out
    return out.getvalue()

def get_input_matrices(fs):
    """
    @return: M, R
    """
    # get the positive strict lower triangular part of the S matrix
    L = []
    for i, line in enumerate(fs.lowtri):
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
    # check reducibility of the exchangeability
    if not MatrixUtil.is_symmetric_irreducible(S):
        raise ValueError('exchangeability is not irreducible')
    # get the mutation rate matrix
    M = S * mut_distn * fs.mutscale
    M -= np.diag(np.sum(M, axis=1))
    # check sign symmetry and irreducibility
    if not MatrixUtil.is_symmetric_irreducible(np.sign(M)):
        msg = 'mutation rate matrix is not sign symmetric irreducible'
        raise ValueError(msg)
    # get the mutation selection balance rate matrix
    R = mrate.to_gtr_halpern_bruno(M, mutsel_distn)
    # check sign symmetry and irreducibility
    if not MatrixUtil.is_symmetric_irreducible(np.sign(R)):
        msg = 'mut-sel balance rate matrix is not sign symmetric irreducible'
        raise ValueError(msg)
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

