"""
This module is about general finite-state continuous-time Markov processes.

The python packages numpy and scipy are used throughout.
At some point I should probably make separate modules
to emphasize the distinction between reversible and general
continuous-time Markov processes.
"""

import math
import unittest
import itertools

import numpy as np
import scipy
from scipy import linalg

import bernoulli
import graph
import MatrixUtil
from MatrixUtil import ndot
import StatsUtil

def sample_distn(n):
    v = np.random.exponential(1, n)
    return v / np.sum(v)

def get_commute_distance_matrix(R, v):
    """
    @param R: reversible rate matrix
    """
    n = len(R)
    psi = np.sqrt(v)
    R_sim = (R.T * psi).T / psi
    R_sim_pinv = scipy.linalg.pinv(R_sim)
    myouter = np.outer(np.ones(n), np.diag(R_sim_pinv))
    D = 2 * R_sim_pinv - myouter - myouter.T
    return D

def get_path_rate_matrix(nstates):
    """
    This is a 3-state path rate matrix.
    The stationary distribution is uniform.
    The sparsity structure is also known as a path graph.
    The matrix is normalized by expected rate.
    """
    A = np.zeros((nstates, nstates))
    for i in range(nstates-1):
        A[i,i+1] = 1.0
        A[i+1,i] = 1.0
    Q = -A
    for i in range(nstates):
        Q[i, i] = -np.sum(Q[i])
    return Q / Q_to_expected_rate(Q)

def get_dense_sequence_rate_matrix(nresidues, nsites):
    """
    Create an reversible rate matrix with uniform stationary distribution.
    Each sequences changes to each other sequence at the same rate.
    The sparsity structure is also known as a complete graph.
    The matrix is normalized by expected rate.
    @param nresidues: for example 4 for DNA or 20 for amino acids
    @param nsites: jointly consider this many sites
    """
    nstates = nresidues**nsites
    R = np.ones((nstates, nstates))
    for i in range(nstates):
        R[i, i] = -(nstates - 1)
    uniform_pi = np.reciprocal(nstates * np.ones(nstates))
    expected_rate = -sum(uniform_pi[i] * R[i, i] for i in range(nstates))
    return R / expected_rate

def get_sparse_sequence_rate_matrix(nresidues, nsites):
    """
    Create an reversible rate matrix with uniform stationary distribution.
    Sites change independently.
    The sparsity structure is also known as a Hamming graph.
    The matrix is normalized by expected rate.
    @param nresidues: for example 4 for DNA or 20 for amino acids
    @param nsites: jointly consider this many sites
    """
    nstates = nresidues**nsites
    R = np.zeros((nstates, nstates))
    for alpha in itertools.product(range(nresidues), repeat=nsites):
        for beta in itertools.product(range(nresidues), repeat=nsites):
            alpha_index = sum(alpha[i]*(nresidues ** i) for i in range(nsites))
            beta_index = sum(beta[i]*(nresidues ** i) for i in range(nsites))
            hamming_dist = sum(1 for a, b in zip(alpha, beta) if a != b)
            if hamming_dist == 1:
                R[alpha_index, beta_index] = 1
    for i in range(nstates):
        R[i, i] = -np.sum(R[i])
    uniform_pi = np.reciprocal(nstates * np.ones(nstates))
    expected_rate = -sum(uniform_pi[i] * R[i, i] for i in range(nstates))
    return R / expected_rate


########################################################################
## These are rate matrix transformations which preserve detailed balance.

def to_gtr_balanced_known_distn(M, p, v):
    """
    @param M: time-reversible rate matrix
    @param p: current stationary distribution
    @param v: target stationary distribution
    @return: a time-reversible rate matrix
    """
    n = len(v)
    R = M.copy()
    # adjust the entries of the rate matrix
    for a in range(n):
        for b in range(n):
            tau = (v[b] / p[b]) / (v[a] / p[a])
            R[a, b] *= math.sqrt(tau)
    # reset the diagonal entries of the rate matrix
    R -= np.diag(np.sum(R, axis=1))
    return R

def to_gtr_balanced(M, v):
    """
    @param M: time-reversible rate matrix
    @param v: target stationary distribution
    @return: a time-reversible rate matrix
    """
    p = R_to_distn(M)
    return to_gtr_balanced_known_distn

def to_gtr_hb_known_energies(M, u, v):
    """
    The hb in the function name stands for Halpern-Bruno.
    @param M: a time-reversible rate matrix
    @param u: energies of states of M
    @param v: target energies
    @return: a time-reversible rate matrix
    """
    n = len(v)
    R = M.copy()
    # adjust the entries of the rate matrix
    for a in range(n):
        for b in range(n):
            if a != b:
                du = u[b] - u[a]
                dv = v[b] - v[a]
                log_tau = du - dv
                coeff = bernoulli.bgf(-log_tau)
                R[a, b] *= coeff
    # reset the diagonal entries of the rate matrix
    R -= np.diag(np.sum(R, axis=1))
    return R


def to_gtr_halpern_bruno_known_distn(M, p, v):
    """
    @param M: a time-reversible rate matrix
    @param p: current stationary distribution
    @param v: target stationary distribution
    @return: a time-reversible rate matrix
    """
    n = len(v)
    R = M.copy()
    # adjust the entries of the rate matrix
    for a in range(n):
        for b in range(n):
            if a != b:
                if p[a] == p[b]:
                    tau = v[b] / v[a]
                elif v[a] == v[b]:
                    tau = p[a] / p[b]
                else:
                    tau = (v[b] / p[b]) / (v[a] / p[a])
                if np.isnan(tau):
                    print a, b, v[b], p[b], v[a], p[a]
                    raise ValueError('tau is nan')
                if not tau:
                    coeff = 0
                elif np.allclose(tau, 1):
                    coeff = 1
                else:
                    coeff = math.log(tau) / (1 - 1/tau)
                R[a, b] *= coeff
    # reset the diagonal entries of the rate matrix
    R -= np.diag(np.sum(R, axis=1))
    return R

def to_gtr_halpern_bruno(M, v):
    """
    @param M: a time-reversible rate matrix
    @param v: target stationary distribution
    @return: a time-reversible rate matrix
    """
    p = R_to_distn(M)
    return to_gtr_halpern_bruno_known_distn(M, p, v)


###########################
## These are misc functions.

def expm_spectral(R, t):
    """
    This is for testing expm_diff_spectral only.
    You should use scipy.linalg.expm instead.
    """
    n = len(R)
    v = R_to_distn(R)
    S = symmetrized(R)
    w, U = np.linalg.eigh(S)
    P = np.zeros_like(R)
    for i in range(n):
        for j in range(n):
            for k in range(n):
                a = (v[j] / v[i])**0.5
                b = U[i, k] * U[j, k]
                c = math.exp(t * w[k])
                P[i, j] += a * b * c
    return P

def expm_diff_spectral(R, t):
    """
    Get the rates of change of transition probabilities at time t.
    Use the spectral representation.
    @return: entrywise derivative of transition matrix at time t
    """
    n = len(R)
    v = R_to_distn(R)
    S = symmetrized(R)
    w, U = np.linalg.eigh(S)
    P_diff = np.zeros_like(R)
    for i in range(n):
        for j in range(n):
            for k in range(n):
                a = (v[j] / v[i])**0.5
                b = U[i, k] * U[j, k]
                c = w[k] * math.exp(t * w[k])
                P_diff[i, j] += a * b * c
    return P_diff

def _R_to_eigenpair(R):
    n = len(R)
    Wl, Vl = scipy.linalg.eig(R, left=True, right=False)
    val_vec_pairs = [(abs(Wl[i]), Vl[:,i]) for i in range(n)]
    r_recip, fiedler = sorted(val_vec_pairs)[1]
    return r_recip, fiedler

def R_to_fiedler(R):
    r_recip, fiedler = _R_to_eigenpair(R)
    return fiedler

def R_to_relaxation_time_obsolete(R):
    """
    This fails when the corresponding eigenvalue is indistinct.
    """
    r_recip, fiedler = _R_to_eigenpair(R)
    return 1 / r_recip

def R_to_relaxation_time(R):
    """
    This assumes a reversible irreducible rate matrix.
    """
    # get the abs eigenvalue directly
    W = np.linalg.eigvals(R)
    abs_eigenvalue = sorted(abs(w) for w in W)[1]
    # get the abs eigenvalue using symmetrization
    W_h = np.linalg.eigvalsh(symmetrized(R))
    abs_eigenvalue_h = sorted(abs(w) for w in W_h)[1]
    # check that the absolute values of the eigenvalues is the same
    if not np.allclose(abs_eigenvalue, abs_eigenvalue_h):
        raise ValueError(
                'relaxation time computation error: %f != %f\n'
                'plain spectrum: %s\n'
                'symmetrized spectrum: %s' % (
                    abs_eigenvalue, abs_eigenvalue_h, W, W_h))
    # return the relaxation time
    return 1 / abs_eigenvalue

def R_to_relaxation_time_experimental(R):
    """
    This assumes a reversible irreducible rate matrix.
    """
    # get the abs eigenvalue directly
    W = scipy.linalg.eigvals(R)
    #TODO just use the second entry of W
    abs_eigenvalue = sorted(abs(w) for w in W)[1]
    # get the abs eigenvalue using symmetrization
    W_h = np.linalg.eigvalsh(symmetrized(R))
    abs_eigenvalue_h = sorted(abs(w) for w in W_h)[1]
    # check that the absolute values of the eigenvalues is the same
    if not np.allclose(abs_eigenvalue, abs_eigenvalue_h):
        raise ValueError(
                'relaxation time computation error: %f != %f' % (
                    abs_eigenvalue, abs_eigenvalue_h))
    # return the relaxation time
    return 1 / abs_eigenvalue

def R_to_distn_nonspectral(R):
    """
    The rate matrix must be irreducible and reversible.
    It is not necessarily symmetric.
    If the rate matrix is symmetric then this function is overkill
    because the stationary distribution would be uniform.
    """
    nstates = len(R)
    V = set(range(nstates))
    E = set()
    for i in range(nstates):
        for j in range(i):
            if R[i, j]:
                if not R[j, i]:
                    raise MatrixUtil.MatrixError('the matrix is not reversible')
                edge = frozenset((i,j))
                E.add(edge)
    nd = graph.g_to_nd(V, E)
    # construct an arbitrary rooted spanning tree of the states
    V_component, D_component = graph.nd_to_dag_component(nd, 0)
    if V_component != V:
        raise MatrixUtil.MatrixError('the matrix is not irreducible')
    # compute the stationary probabilities relative to the first state
    weights = [None] * nstates
    v_to_children = graph.dag_to_cd(V_component, D_component)
    preorder_states = graph.topo_sort(V_component, D_component)
    weights[preorder_states[0]] = 1.0
    for parent in preorder_states:
        for child in v_to_children[parent]:
            ratio = R[parent, child] / R[child, parent]
            weights[child] = weights[parent] * ratio
    total = sum(weights)
    return np.array(weights) / total

def P_to_distn(R):
    """
    @param P: transition matrix
    @return: stationary distribution
    """
    n = len(R)
    Wl, Vl = scipy.linalg.eig(R, left=True, right=False)
    val_vec_pairs = [(abs(Wl[i]), Vl[:,i]) for i in range(n)]
    dummy, pi_eigenvector = max(val_vec_pairs)
    total = np.sum(pi_eigenvector)
    pi_arr = np.array([v/total for v in pi_eigenvector])
    return pi_arr

def R_to_distn(R):
    """
    This seems to have a problem when eigenvalues are not distinct.
    The problem is that it tries to compare the eigenvectors.
    @param R: rate matrix
    @return: stationary distribution
    """
    n = len(R)
    Wl, Vl = scipy.linalg.eig(R, left=True, right=False)
    val_vec_pairs = [(abs(Wl[i]), Vl[:,i]) for i in range(n)]
    dummy, pi_eigenvector = min(val_vec_pairs)
    total = np.sum(pi_eigenvector)
    pi_arr = np.array([v/total for v in pi_eigenvector])
    return pi_arr

def Q_to_expected_rate(Q):
    """
    The rate matrix Q should be irreducible.
    But it does not have to be reversible.
    The term "expected rate" is preferred by Jeff Thorne
    and is used in the wikipedia article about
    models of DNA evolution.
    The term is also used in source code in PAL
    by Korbinian Strimmer and Alexei Drummond
    as the "expected number of substitutions".
    This quantitiy is also called "total rate of change per unit time",
    for example in the paper
    "Toward Extracting All Phylogenetic Information
    from Matrices of Evolutionary Distances"
    by Sebastien Roch.
    @param Q: a rate matrix
    """
    n = len(Q)
    distn = R_to_distn(Q)
    expected_rate = 0.0
    for i in range(n):
        expected_rate -= distn[i] * Q[i, i]
    return expected_rate

def symmetrized(R):
    """
    Get the symmetrized matrix of a reversible markov process.
    This returns a symmetric matrix that is not a rate matrix
    because rows do not sum to zero.
    The returned matrix should be similar to R
    in the sense of linear algebra matrix similarity.
    """
    v = R_to_distn(R)
    lam = np.diag(np.sqrt(v))
    rlam = np.diag(np.reciprocal(np.sqrt(v)))
    return ndot(lam, R, rlam)

def symmetrized_known_distn(R, v):
    """
    Get the symmetrized matrix of a reversible markov process.
    This returns a symmetric matrix that is not a rate matrix
    because rows do not sum to zero.
    The returned matrix should be similar to R
    in the sense of linear algebra matrix similarity.
    """
    lam = np.diag(np.sqrt(v))
    rlam = np.diag(np.reciprocal(np.sqrt(v)))
    return ndot(lam, R, rlam)

def _holmes_rubin_2002_kernel_function(mu_k, mu_l, T):
    """
    This could be re-implemented more stably.
    @param mu_k: a rate matrix eigenvalue
    @param mu_l: a rate matrix eigenvalue
    @param T: an amount of time
    """
    #if mu_k == mu_l:
    if np.allclose(mu_k, mu_l):
        return T * math.exp(mu_k * T)
    else:
        return (math.exp(mu_k * T) - math.exp(mu_l * T)) / (mu_k - mu_l)

def _holmes_rubin_2002_kernel(w, T):
    """
    @param w: eigenvalues whose corresponding eigenvectors are columns of U
    @param T: an amount of time
    @return: a symmetric matrix
    """
    return np.array([[
        _holmes_rubin_2002_kernel_function(x, y, T) for x in w] for y in w])

def _holmes_rubin_2002_summation(U, a, b, i, K):
    """
    @param U: an orthonormal matrix
    @param a: integer initial state index
    @param b: integer final state index
    @param i: integer query state index
    @param K: a symmetric matrix with eigenvalue and time information
    """
    """
    total = 0
    for k in range(len(U)):
        for l in range(len(U)):
            total += U[a,k]*U[i,k]*U[i,l]*U[b,l]*K[k,l]
    return total
    """
    u = U[a] * U[i]
    v = U[b] * U[i]
    return ndot(u, K, v)


def get_endpoint_conditioned_expected_occupancy(R, v, a, b, T):
    """
    Holmes and Rubin 2002.
    @param R: rate matrix
    @param v: stationary distribution
    @param a: integer state index of initial state
    @param b: integer state index of final state
    @param T: elapsed time
    @return: endpoint conditioned expected amount of time spent in each state
    """
    n = len(v)
    psi = np.sqrt(v)
    S = (R.T * psi).T / psi
    MatrixUtil.assert_symmetric(S)
    w, U = scipy.linalg.eigh(S)
    if not np.allclose(np.dot(U, U.T), np.eye(n)):
        raise Exception('U should be orthogonal')
    P = scipy.linalg.expm(T*R)
    # the Mab is Holmes and Rubin 2002 notation
    Mab = (psi[b] / psi[a]) * np.sum(U[a] * U[b] * np.exp(T*w))
    if not np.allclose(P[a,b], Mab):
        raise Exception('not close: %s %s' % (P[a,b], Mab))
    coeff = (psi[b] / psi[a]) / Mab
    K = _holmes_rubin_2002_kernel(w, T)
    occupancy = coeff * np.array([
        _holmes_rubin_2002_summation(U, a, b, i, K) for i in range(n)])
    if not np.allclose(T, np.sum(occupancy)):
        raise Exception(
                'the expectected occupancy times should add up '
                'to the total time')
    return occupancy

def get_hobolth_eceo(R, v, a, b, T, nmax):
    """
    The eceo means endpoint cnoditioned expected occupancy.
    Most of the function arguments are the same as those of the more 
    verbosely named function.
    @param nmax: truncation of an infinite summation
    """
    accum = np.zeros(len(v))
    mu = np.max(-np.diag(R))
    X = np.eye(len(v)) + R / mu
    for n in range(nmax+1):
        coeff = (T / (n+1)) * math.exp(StatsUtil.poisson_log_pmf(n, mu*T))
        #print 'coeff:', coeff
        for alpha in range(len(v)):
            conditional_sum = 0
            for i in range(n+1):
                prefix = np.linalg.matrix_power(X, i)[a, alpha]
                suffix = np.linalg.matrix_power(X, n-i)[alpha, b]
                conditional_sum += prefix * suffix
                #print 'conditional sum:', conditional_sum
            accum[alpha] += coeff * conditional_sum
    return accum / scipy.linalg.expm(R*T)[a,b]


class TestMrate(unittest.TestCase):

    def test_expm(self):
        M = np.array([
            [-1.55273124,  0.40323905,  0.90129456,  0.24819763],
            [ 2.41191856, -5.14753325,  2.29550348,  0.44011122],
            [ 0.82711917,  0.3521918,  -1.44057109,  0.26126012],
            [ 3.01693453,  0.89439764,  3.46050956, -7.37184173]])
        t = 0.3
        observed = scipy.linalg.expm(M * t)
        expected = expm_spectral(M, t)
        self.assertTrue(np.allclose(observed, expected))

    def test_sample_distn(self):
        n = 5
        v = sample_distn(n)
        self.assertEqual(v.shape, (n,))
        self.assertTrue(np.allclose(np.sum(v), 1))
        self.assertTrue(np.all(v > 0))

    def test_holmes_rubin(self):
        v = np.array([.4, .2, .2, .2])
        R = np.array([
            [-3, 1, 1, 1],
            [2, -4, 1, 1],
            [2, 1, -4, 1],
            [2, 1, 1, -4]], dtype=float)
        a = 0
        b = 1
        T = 1.0
        spectral = get_endpoint_conditioned_expected_occupancy(R, v, a, b, T)
        nmax = 60
        uniformized = get_hobolth_eceo(R, v, a, b, T, nmax)
        self.assertTrue(np.allclose(spectral, uniformized))


if __name__ == '__main__':
    unittest.main()
