r"""
How exactly does Sella-Hirsh fitness reduction affect the connectedness.

We consider very restricted mutation models
whose nonzero rates are all the same,
and whose sparsity structures are like lattices or
complete graphs, or hypercubes, or hamming graphs.
Then we reduce the fitness of a single state towards zero
using the Sella-Hirsh approximation and look at the effect
on the mutation-selection balance rate matrix.
In particular we want to know how the expected rate
and the relaxation time are affected
as the stationary probability of the unfit state approaches zero.
Is the selective removal of this state more like
vertex deletion, or is it more like removal by Schur complementation?
The simple selection option is like the
Knudsen-Miyamoto selection or in other words
the Goldman and Whelan selection with \( f = \frac{1}{2} \).
The complicated selection option is the one with the logarithm
that Jeff Thorne likes to use and which was also
used by Sella and Hirsh, and by Halpern and Bruno.
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

def get_form():
    """
    @return: the body of a form
    """
    # define the form objects
    form_objects = [
            Form.Float('p_mid',
                'mutation-selection probability for the down-weighted state',
                '0.01', low_exclusive=0, high_exclusive=1),
            Form.RadioGroup('mut_type', 'mutation graph structure', [
                Form.RadioItem('mut_path_2', 'path (2 states)'),
                Form.RadioItem('mut_path_3_mid',
                    'path (3 states; middle is down-weighted)'),
                Form.RadioItem('mut_path_3_end',
                    'path (3 states; end is down-weighted)'),
                Form.RadioItem('mut_complete_3', 'complete graph (3 states)'),
                Form.RadioItem('mut_complete_4', 'complete graph (4 states)'),
                Form.RadioItem('mut_hyper_2_2', 'square H(2, 2)', True),
                Form.RadioItem('mut_hyper_2_3', 'cube H(2, 3)')]),
            Form.RadioGroup('sel_type', 'selection approximation', [
                Form.RadioItem('sel_sqrt', 'simple'),
                Form.RadioItem('sel_log', 'sophisticated', True)])]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_detailed_balance_error(Q):
    """
    @param Q: a rate matrix
    @return: a number that should be near zero if detailed balance is satisfied
    """
    p = mrate.R_to_distn(Q)
    errors = []
    nstates = len(Q)
    for i in range(nstates):
        for j in range(nstates):
            error = p[i] * Q[i, j] - p[j] * Q[j, i]
            errors.append(error)
    return min(abs(x) for x in errors)

def get_rate_matrix_summary(Q):
    out = StringIO()
    Q_v = mrate.R_to_distn(Q)
    Q_r = mrate.R_to_total_rate(Q)
    Q_t = mrate.R_to_relaxation_time(Q)
    print >> out, 'rate matrix:'
    print >> out, Q
    print >> out
    print >> out, 'this should be near zero for detailed balance:'
    print >> out, get_detailed_balance_error(Q)
    print >> out
    print >> out, 'computed stationary distribution:'
    print >> out, Q_v
    print >> out
    print >> out, 'expected rate:'
    print >> out, Q_r
    print >> out
    print >> out, 'relaxation time'
    print >> out, Q_t
    print >> out
    print >> out, '(expected rate) * (relaxation time):'
    print >> out, Q_r * Q_t
    print >> out
    print >> out
    return out.getvalue().rstrip()


def do_mut_hyper_2_3(fs, to_gtr):
    out = StringIO()
    # define the path mutation rate matrix
    M = mrate.get_sparse_sequence_rate_matrix(2, 3)
    print >> out, '*** mutation rate matrix (8-state cube) ***'
    print >> out
    print >> out, get_rate_matrix_summary(M)
    print >> out
    print >> out
    # kill the last state by natural selection
    p_other = (1 - fs.p_mid)/7
    p_target = [p_other]*7 + [fs.p_mid]
    Q = to_gtr(M, p_target)
    print >> out, '*** mutation-selection balance ***'
    print >> out
    print >> out, get_rate_matrix_summary(Q)
    print >> out
    print >> out
    # define a reference mutation rate matrix
    R = mrate.get_sparse_sequence_rate_matrix(2, 3)
    nstates = 7
    M = np.zeros((nstates, nstates))
    for i in range(nstates):
        for j in range(nstates):
            if i != j:
                M[i, j] = R[i, j]
    M -= np.diag(np.sum(M, axis=1))
    M /= mrate.R_to_total_rate(M)
    print >> out, '*** reference mutation rate matrix (corner removed) ***'
    print >> out
    print >> out, get_rate_matrix_summary(M)
    print >> out
    print >> out
    return out.getvalue().rstrip()

def do_mut_hyper_2_2(fs, to_gtr):
    out = StringIO()
    # define the path mutation rate matrix
    M = mrate.get_sparse_sequence_rate_matrix(2, 2)
    print >> out, '*** mutation rate matrix (4-state square) ***'
    print >> out
    print >> out, get_rate_matrix_summary(M)
    print >> out
    print >> out
    # kill the last state by natural selection
    p_other = (1 - fs.p_mid)/3
    p_target = (p_other, p_other, p_other, fs.p_mid)
    Q = to_gtr(M, p_target)
    print >> out, '*** mutation-selection balance ***'
    print >> out
    print >> out, get_rate_matrix_summary(Q)
    print >> out
    print >> out
    # define a reference mutation rate matrix
    M = mrate.get_path_rate_matrix(3)
    print >> out, '*** reference mutation rate matrix (3-state path) ***'
    print >> out
    print >> out, get_rate_matrix_summary(M)
    print >> out
    print >> out
    return out.getvalue().rstrip()

def do_mut_path_2(fs, to_gtr):
    out = StringIO()
    # define the path mutation rate matrix
    M = mrate.get_path_rate_matrix(2)
    print >> out, '*** mutation rate matrix (2-state path) ***'
    print >> out
    print >> out, get_rate_matrix_summary(M)
    print >> out
    print >> out
    # kill the last state by natural selection
    p_target = ((1 - fs.p_mid), fs.p_mid)
    Q = to_gtr(M, p_target)
    print >> out, '*** mutation-selection balance ***'
    print >> out
    print >> out, get_rate_matrix_summary(Q)
    print >> out
    print >> out
    return out.getvalue().rstrip()

def do_mut_path_3_mid(fs, to_gtr):
    out = StringIO()
    # define the path mutation rate matrix
    M = mrate.get_path_rate_matrix(3)
    print >> out, '*** mutation rate matrix (3-state path) ***'
    print >> out
    print >> out, get_rate_matrix_summary(M)
    print >> out
    print >> out
    # kill the middle state by natural selection
    p_target = ((1 - fs.p_mid)/2, fs.p_mid, (1 - fs.p_mid)/2)
    Q = to_gtr(M, p_target)
    print >> out, '*** mutation-selection balance ***'
    print >> out
    print >> out, get_rate_matrix_summary(Q)
    print >> out
    print >> out
    # define a reference mutation rate matrix that is a scaled schur complement
    M = mrate.get_path_rate_matrix(2)
    print >> out, '*** reference mutation rate matrix (2-state path) ***'
    print >> out
    print >> out, get_rate_matrix_summary(M)
    print >> out
    print >> out
    return out.getvalue().rstrip()

def do_mut_path_3_end(fs, to_gtr):
    out = StringIO()
    # define the path mutation rate matrix
    M = mrate.get_path_rate_matrix(3)
    print >> out, '*** mutation rate matrix (3-state path) ***'
    print >> out
    print >> out, get_rate_matrix_summary(M)
    print >> out
    print >> out
    # kill the middle state by natural selection
    p_target = ((1 - fs.p_mid)/2, (1 - fs.p_mid)/2, fs.p_mid)
    Q = to_gtr(M, p_target)
    print >> out, '*** mutation-selection balance ***'
    print >> out
    print >> out, get_rate_matrix_summary(Q)
    print >> out
    print >> out
    # define a reference mutation rate matrix that is a scaled schur complement
    M = mrate.get_path_rate_matrix(2)
    print >> out, '*** reference mutation rate matrix (2-state path) ***'
    print >> out
    print >> out, get_rate_matrix_summary(M)
    print >> out
    print >> out
    return out.getvalue().rstrip()

def do_mut_complete_3(fs, to_gtr):
    out = StringIO()
    # define the path mutation rate matrix
    M = mrate.get_dense_sequence_rate_matrix(3, 1)
    print >> out, '*** mutation rate matrix (3-state complete) ***'
    print >> out
    print >> out, get_rate_matrix_summary(M)
    print >> out
    print >> out
    # kill the middle state by natural selection
    p_other = (1 - fs.p_mid)/2
    p_target = (p_other, p_other, fs.p_mid)
    Q = to_gtr(M, p_target)
    print >> out, '*** mutation-selection balance ***'
    print >> out
    print >> out, get_rate_matrix_summary(Q)
    print >> out
    print >> out
    # define a reference mutation rate matrix
    M = mrate.get_path_rate_matrix(2)
    print >> out, '*** reference mutation rate matrix (2-state path) ***'
    print >> out
    print >> out, get_rate_matrix_summary(M)
    print >> out
    print >> out
    return out.getvalue().rstrip()

def do_mut_complete_4(fs, to_gtr):
    out = StringIO()
    # define the path mutation rate matrix
    M = mrate.get_dense_sequence_rate_matrix(4, 1)
    print >> out, '*** mutation rate matrix (4-state complete) ***'
    print >> out
    print >> out, get_rate_matrix_summary(M)
    print >> out
    print >> out
    # kill the middle state by natural selection
    p_other = (1 - fs.p_mid)/3
    p_target = (p_other, p_other, p_other, fs.p_mid)
    Q = to_gtr(M, p_target)
    print >> out, '*** mutation-selection balance ***'
    print >> out
    print >> out, get_rate_matrix_summary(Q)
    print >> out
    print >> out
    # define a reference mutation rate matrix that is a scaled schur complement
    M = mrate.get_dense_sequence_rate_matrix(3, 1)
    print >> out, '*** reference mutation rate matrix (3-state complete) ***'
    print >> out
    print >> out, get_rate_matrix_summary(M)
    print >> out
    print >> out
    return out.getvalue().rstrip()


def get_response_content(fs):
    out = StringIO()
    np.set_printoptions(linewidth=200)
    # define the selection type
    if fs.sel_sqrt:
        to_gtr = mrate.to_gtr_balanced
    elif fs.sel_log:
        to_gtr = mrate.to_gtr_halpern_bruno
    # decide which mutation process to analyze
    if fs.mut_path_2:
        print >> out, do_mut_path_2(fs, to_gtr)
    elif fs.mut_path_3_mid:
        print >> out, do_mut_path_3_mid(fs, to_gtr)
    elif fs.mut_path_3_end:
        print >> out, do_mut_path_3_end(fs, to_gtr)
    elif fs.mut_complete_3:
        print >> out, do_mut_complete_3(fs, to_gtr)
    elif fs.mut_complete_4:
        print >> out, do_mut_complete_4(fs, to_gtr)
    elif fs.mut_hyper_2_2:
        print >> out, do_mut_hyper_2_2(fs, to_gtr)
    elif fs.mut_hyper_2_3:
        print >> out, do_mut_hyper_2_3(fs, to_gtr)
    return out.getvalue()

