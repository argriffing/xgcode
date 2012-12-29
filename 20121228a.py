"""
Check eigendecompositions of block triangular rate matrices.
"""

from StringIO import StringIO

import numpy as np
import scipy.linalg

import Form
import FormOut
from MatrixUtil import ndot

def get_form():
    """
    @return: the body of a form
    """
    return [
            Form.RadioGroup('mattype', 'predefined rate matrix', [
                Form.RadioItem('plain', 'plain 4x4 integer rate matrix', True),
                Form.RadioItem('rand', 'random rate matrix'),
                ]),
            ]

def get_form_out():
    return FormOut.Report()

def assert_reversibility(Q, v):
    if Q.ndim != 2:
        raise Exception
    if v.ndim != 1:
        raise Exception
    nstates = v.shape[0]
    if Q.shape != (nstates, nstates):
        raise Exception
    if not np.all(np.greater_equal(v, 0)):
        raise Exception
    if not np.allclose(np.sum(v), 1):
        raise Exception
    if not np.all(np.less(np.diag(Q), 0)):
        raise Exception
    if not np.all(np.greater_equal(Q - np.diag(np.diag(Q)), 0)):
        raise Exception
    if not np.allclose(np.dot(v, Q), 0):
        raise Exception
    for i in range(nstates):
        for j in range(nstates):
            if not np.allclose(v[i] * Q[i, j], v[j] * Q[j, i]):
                raise Exception

def sample_distn(nstates):
    v = np.exp(np.random.randn(nstates))
    return v / np.sum(v)

def get_reversible_rate_matrix(nstates):
    """
    Return a random reversible rate matrix and its stationary distribution.
    This function is used to help construct structured rate matrices
    with transient states.
    @param nstates: number of states
    @return: Q, v
    """
    M = np.exp(np.random.randn(nstates, nstates))
    pre_Q = 0.5 * (M + M.T)
    v = sample_distn(nstates)
    v_sqrt = np.sqrt(v)
    pre_Q = (pre_Q.T / v_sqrt).T * v_sqrt
    Q = pre_Q - np.diag(np.sum(pre_Q, axis=1))
    assert_reversibility(Q, v)
    return Q, v

def get_plain_rate_matrix():
    pre_Q = np.array([
        [0, 1, 1, 0],
        [1, 0, 0, 1],
        [0, 0, 0, 1],
        [0, 0, 1, 0],
        ], dtype=float)
    Q = pre_Q - np.diag(np.sum(pre_Q, axis=1))
    v_trans = 0.5 * np.ones(2)
    v_recur = 0.5 * np.ones(2)
    return Q, v_trans, v_recur

def get_random_structured_rate_matrix():
    """
    The first few states are transient.
    """
    nstates_block = 3
    Q_trans, v_trans = get_reversible_rate_matrix(nstates_block)
    Q_recur, v_recur = get_reversible_rate_matrix(nstates_block)
    outflow = np.exp(np.random.randn(nstates_block))
    Q11 = Q_trans
    Q12 = np.diag(outflow)
    Q21 = np.zeros((nstates_block, nstates_block))
    Q22 = Q_recur
    pre_Q = np.array(np.bmat([[Q11, Q12], [Q21, Q22]]))
    Q = pre_Q - np.diag(np.sum(pre_Q, axis=1))
    v = np.hstack((np.zeros_like(v_trans), v_recur))
    assert_reversibility(Q, v)
    return Q, v_trans, v_recur

def get_response_content(fs):
    if fs.plain:
        Q, v_trans, v_recur = get_plain_rate_matrix()
    elif fs.rand:
        Q, v_trans, v_recur = get_random_structured_rate_matrix()
    else:
        raise Exception
    nstates = Q.shape[0]
    w, vl, vr = scipy.linalg.eig(Q, left=True, right=True)
    vl_inv = scipy.linalg.inv(vl)
    vr_inv = scipy.linalg.inv(vr)
    #
    # do weird things with the sylvester equation
    n = nstates / 2
    A_syl = Q[:n, :n]
    B_syl = -Q[n:, n:]
    Q_syl = -Q[:n, n:]
    X = scipy.linalg.solve_sylvester(A_syl, B_syl, Q_syl)
    T = np.array(np.bmat([
        [np.eye(n), X],
        [np.zeros((n, n)), np.eye(n)],
        ]))
    T_inv = scipy.linalg.inv(T)
    #
    # do stuff with the stationary distributions of the separate processes
    v_trans_recur = np.hstack((v_trans, v_recur))
    D_sqrt = np.diag(np.sqrt(v_trans_recur))
    D_sqrt_recip = np.diag(np.reciprocal(np.sqrt(v_trans_recur)))
    block_diag = ndot(D_sqrt, T_inv, Q, T, D_sqrt_recip)
    w, U = scipy.linalg.eigh(block_diag)
    w_full = np.diag(ndot(U.T, D_sqrt, T_inv, Q, T, D_sqrt_recip, U))
    #
    np.set_printoptions(
            linewidth=1000000,
            threshold=1000000,
            )
    out = StringIO()
    print >> out, 'Q:'
    print >> out, Q
    print >> out
    print >> out, 'w:'
    print >> out, w
    print >> out
    print >> out, 'vl:'
    print >> out, vl
    print >> out
    print >> out, 'vl.T:'
    print >> out, vl.T
    print >> out
    print >> out, 'inv(vl):'
    print >> out, vl_inv
    print >> out
    print >> out, 'vr:'
    print >> out, vr
    print >> out
    print >> out, 'vr.T:'
    print >> out, vr.T
    print >> out
    print >> out, 'inv(vr):'
    print >> out, vr_inv
    print >> out
    print >> out, 'inv(vl).T w vl.T:'
    print >> out, np.dot(vl_inv.T, np.dot(np.diag(w), vl.T))
    print >> out
    print >> out, 'vr w inv(vr):'
    print >> out, np.dot(vr, np.dot(np.diag(w), vr_inv))
    print >> out
    print >> out
    print >> out, 'sylvester equation stuff...'
    print >> out
    print >> out, 'X:'
    print >> out, X
    print >> out
    print >> out, 'T:'
    print >> out, T
    print >> out
    print >> out, 'inv(T):'
    print >> out, T_inv
    print >> out
    print >> out, 'inv(T) Q T:'
    print >> out, np.dot(T_inv, np.dot(Q, T))
    print >> out
    print >> out, 'U.T D^(1/2) inv(T) Q T D^(-1/2) U:'
    print >> out, ndot(U.T, D_sqrt, T_inv, Q, T, D_sqrt_recip, U)
    print >> out
    print >> out, 'expm(Q):'
    print >> out, scipy.linalg.expm(Q)
    print >> out
    print >> out, 'T D^-1/2 U exp(w) U.T D^1/2 T^-1'
    print >> out, ndot(
            T, D_sqrt_recip, U, np.diag(np.exp(w_full)), U.T, D_sqrt, T_inv)
    print >> out
    return out.getvalue()

