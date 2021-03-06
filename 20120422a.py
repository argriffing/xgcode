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
import cheeger
import MatrixUtil
from MatrixUtil import ndot

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
                    'infotype' : 'info_mut'}),
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
    return FormOut.Html('summary')

class RateMatrixSummary:
    def __init__(self, Q):
        """
        @param Q: rate matrix
        """
        # define intermediate variables
        v = mrate.R_to_distn(Q)
        n = len(v)
        psi = np.sqrt(v)
        c_low, c_mid, c_high = cheeger.get_cheeger_bounds(Q, v)
        # define member variables to summarize the rate matrix
        self.rate_matrix = Q
        self.exch_matrix = Q / v
        if not np.allclose(self.exch_matrix, self.exch_matrix.T):
            print self.exch_matrix
            raise ValueError('expected symmetry')
        self.sim_sym_matrix = np.outer(psi, 1/psi) * Q
        if not np.allclose(self.sim_sym_matrix, self.sim_sym_matrix.T):
            print self.sim_sym_matrix
            raise ValueError('expected symmetry')
        self.distn = v
        self.distn_shannon_entropy = -ndot(np.log(v), v)
        self.distn_logical_entropy = ndot(v, 1-v)
        self.expected_rate = -ndot(np.diag(Q), v)
        self.spectrum = scipy.linalg.eigvalsh(self.sim_sym_matrix)
        self.spectral_gap = -self.spectrum[-2]
        self.isoperimetric_low = c_low
        self.isoperimetric_constant = c_mid
        self.isoperimetric_high = c_high
        self.trace_bound_high = -sum(np.diag(Q)) / (n-1)

def get_html_table(summaries):
    """
    Return text for an html table comparing rate matrix summaries.
    """
    out = StringIO()
    summary_names = [
            'distn shannon entropy',
            'distn logical entropy',
            'expected rate',
            'spectral gap',
            'isoperimetric lower bound on spectral gap',
            'isoperimetric constant',
            'isoperimetric upper bound on spectral gap',
            'trace upper bound on spectral gap',
            ]
    data_columns = []
    for x in summaries:
        col = [
                x.distn_shannon_entropy,
                x.distn_logical_entropy,
                x.expected_rate,
                x.spectral_gap,
                x.isoperimetric_low,
                x.isoperimetric_constant,
                x.isoperimetric_high,
                x.trace_bound_high,
                ]
        data_columns.append(col)
    headers = ('', 'pure mutation', 'mutation selection balance')
    rows = [headers] + zip(*([summary_names] + data_columns))
    print >> out, '<table border="1">'
    for row in rows:
        print >> out, '<tr>'
        for item in row:
            print >> out, '<td>'
            print >> out, item
            print >> out, '</td>'
        print >> out, '</tr>'
    print >> out, '</table>'
    return out.getvalue().rstrip()

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

def get_mi_asymptotics(M, R):
    """
    Return a multiline string with some asymptotics.
    @param M: pure mutation rate matrix
    @param R: mutation-selection balance rate matrix
    @return: multiline string
    """
    out = StringIO()
    # get the stationary distributions
    M_v = mrate.R_to_distn(M)
    R_v = mrate.R_to_distn(R)
    # The shannon entropy of the stationary distribution of the process
    # determines the mutual information at small times.
    M_shannon_entropy = -np.dot(np.log(M_v), M_v)
    R_shannon_entropy = -np.dot(np.log(R_v), R_v)
    if not np.allclose(M_shannon_entropy, R_shannon_entropy):
        print >> out, 'At small enough times'
        if R_shannon_entropy < M_shannon_entropy:
            print >> out, '* pure mutation',
        else:
            print >> out, '* mutation-selection balance',
        print >> out, 'will be more informative'
        print >> out, 'because its stationary distribution has greater',
        print >> out, 'Shannon entropy.'
    else:
        print >> out, 'There is not enough difference between the'
        print >> out, 'Shannon entropies of the stationary distributions'
        print >> out, 'to determine which process'
        print >> out, 'is more informative at times near zero'
    print >> out
    # The spectral gap of the process
    # determines the mutual information at large times.
    M_spectral_gap = sorted(abs(w) for w in scipy.linalg.eigvals(M))[1]
    R_spectral_gap = sorted(abs(w) for w in scipy.linalg.eigvals(R))[1]
    M_cheeg_low, M_cheeg_mid, M_cheeg_high = cheeger.get_cheeger_bounds(M, M_v)
    R_cheeg_low, R_cheeg_mid, R_cheeg_high = cheeger.get_cheeger_bounds(R, R_v)
    if not np.allclose(M_spectral_gap, R_spectral_gap):
        print >> out, 'At large enough times'
        if R_spectral_gap < M_spectral_gap:
            print >> out, '* mutation-selection balance',
        else:
            print >> out, '* pure mutation',
        print >> out, 'will be more informative'
        print >> out, 'because it has a smaller spectral gap.'
        if (R_cheeg_high < M_cheeg_low) or (M_cheeg_high < R_cheeg_low):
            print >> out, 'And also because of the isoperimetric bounds.'
    else:
        print >> out, 'There is not enough difference between the'
        print >> out, 'spectral gaps to determine which process'
        print >> out, 'is more informative at times near infinity'
    print >> out
    # return the text
    return out.getvalue().strip()

def get_heuristics(M, R):
    """
    Return a multiline string with some heuristics.
    The heuristics are independendent of time and of the information variant.
    Greater stationary distribution shannon entropy suggests less saturation.
    Greater stationary distribution logical entropy suggests less saturation.
    Greater expected rate suggests more saturation.
    Greater spectral rate suggests more saturation.
    @param M: pure mutation rate matrix
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
    print >> out, '<html>'
    print >> out, '<body>'
    print >> out
    print >> out, '<pre>'
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
    if fs.info_mut:
        print >> out, 'Mutual information properties',
        print >> out, 'at very small and very large times:'
        print >> out
        print >> out, get_mi_asymptotics(M, R)
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
    print >> out, '</pre>'
    #
    # create the summaries
    summaries = (RateMatrixSummary(M), RateMatrixSummary(R))
    print >> out, get_html_table(summaries)
    print >> out
    print >> out, '<html>'
    print >> out, '<body>'
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
            raise ValueError(
                    'expected %d values on line "%s"' % (i+1, line))
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
        raise ValueError('the inputs do not agree on the state space size')
    # check for sufficient number of states
    if nstates < 2:
        raise ValueError('at least two states are required')
    # check reducibility of the exchangeability
    if not MatrixUtil.is_symmetric_irreducible(S):
        raise ValueError('exchangeability is not irreducible')
    # get the mutation rate matrix
    M = S * mut_distn * fs.mutscale
    M -= np.diag(np.sum(M, axis=1))
    # check sign symmetry and irreducibility
    if not MatrixUtil.is_symmetric_irreducible(np.sign(M)):
        raise ValueError(
                'mutation rate matrix is not sign symmetric irreducible')
    # get the mutation selection balance rate matrix
    R = mrate.to_gtr_halpern_bruno(M, mutsel_distn)
    # check sign symmetry and irreducibility
    if not MatrixUtil.is_symmetric_irreducible(np.sign(R)):
        raise ValueError(
                'mut-sel balance rate matrix '
                'is not sign symmetric irreducible')
    # check the stationary distributions
    mut_distn_observed = mrate.R_to_distn(M)
    if not np.allclose(mut_distn_observed, mut_distn):
        raise ValueError(
                'internal mut stationary distribution computation error')
    mutsel_distn_observed = mrate.R_to_distn(R)
    if not np.allclose(mutsel_distn_observed, mutsel_distn):
        raise ValueError(
                'internal mut-sel stationary distribution computation error')
    # return the values
    return M, R

