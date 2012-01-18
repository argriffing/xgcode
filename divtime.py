import unittest
import math

import numpy as np
import scipy
from scipy import linalg

import mrate
from MatrixUtil import ndot

def sample_distribution(n):
    """
    @param n: number of states
    """
    # Get a nonnegative vector.
    v = np.random.rand(n)
    # Divide the vector by its nonnegative sum.
    distn = v / np.sum(v)
    return distn

def sample_symmetric_rate_matrix(n):
    """
    @param n: number of states
    """
    # Get a nonnegative asymmetric matrix.
    M = np.random.rand(n, n)
    # Symmetrize by adding to the transpose.
    S = M + M.T
    # Subtract row sum from diagonal.
    R = S - np.diag(np.sum(S, axis=1))
    return R

def to_gtr_halpern_bruno(S, v):
    """
    @param S: symmetric rate matrix
    @param v: target stationary distribution
    @return: a rate matrix with the target stationary distribution
    """
    p = mrate.R_to_distn(S)
    # get the number of states
    n = len(v)
    # copy the symmetric rate matrix
    R = S.copy()
    # adjust the entries of the rate matrix
    for a in range(n):
        for b in range(n):
            if a != b:
                # This equation is unnecessarily verbose
                # due to symmetry of S.
                # It should also work for asymmetric input rate matrices.
                tau = (v[b] / p[b]) / (v[a] / p[a])
                if not np.allclose(tau, 1):
                    R[a, b] *= math.log(tau) / (1 - 1/tau)
    # reset the diagonal entries of the rate matrix
    R -= np.diag(np.sum(R, axis=1))
    return R

def get_mutual_information_stable(R, t):
    """
    This is a more stable function.
    @return: unscaled_result, log_of_scaling_factor
    """
    #FIXME under construction
    n = len(R)
    v = mrate.R_to_distn(R)
    S = mrate.symmetrized(R)
    w, U = np.linalg.eigh(S)
    P = np.zeros_like(R)
    accum = 0
    for i in range(n):
        for j in range(n):
            for k in range(n):
                a = (v[j] / v[i])**0.5
                b = U[i, k] * U[j, k]
                c = math.exp(t * w[k])
                P[i, j] += a * b * c
    # compute the unscaled part of log(X(i,j)/(X(i)*X(j)))
    for i in range(n):
        for j in range(n):
            if v[i] and P[i, j]:
                coeff = v[i] * P[i, j]
                numerator = P[i, j]
                denominator = v[j]
                # the problem is that the following log is nearly zero
                value = coeff * math.log(numerator / denominator)
                accum += np.real(value)
    return accum

def get_mutual_information_approx(R, t):
    """
    This is an approximation for large times.
    It can be rewritten using orthogonality.
    """
    n = len(R)
    v = mrate.R_to_distn(R)
    S = mrate.symmetrized(R)
    w, U = np.linalg.eigh(S)
    accum = 0
    for i in range(n):
        for j in range(n):
            b = 0
            for k in range(n-1):
                b += U[i,k]*U[j,k]*math.exp(t*w[k])
            accum += (b * b) / 2
    return accum

def get_mutual_information_approx_b(R, t):
    """
    This is an approximation for large times.
    It has been rewritten using orthogonality.
    """
    n = len(R)
    v = mrate.R_to_distn(R)
    S = mrate.symmetrized(R)
    w, U = np.linalg.eigh(S)
    accum = 0
    for i in range(n):
        for j in range(n):
            for k in range(n-1):
                accum += ((U[i,k]*U[j,k])**2) * math.exp(2*t*w[k]) / 2
    return accum

def get_mutual_information_diff_approx(R, t):
    """
    This is an approximation for large times.
    It can be rewritten using orthogonality.
    """
    n = len(R)
    v = mrate.R_to_distn(R)
    S = mrate.symmetrized(R)
    w, U = np.linalg.eigh(S)
    accum = 0
    for i in range(n):
        for j in range(n):
            b = 0
            for k in range(n-1):
                b += U[i,k]*U[j,k]*math.exp(t*w[k])
            c = 0
            for k in range(n-1):
                c += U[i,k]*U[j,k]*w[k]*math.exp(t*w[k])
            accum += b * c
    return accum

def get_mutual_information_diff_approx_b(R, t):
    """
    This is an approximation for large times.
    It has been rewritten using orthogonality.
    """
    n = len(R)
    v = mrate.R_to_distn(R)
    S = mrate.symmetrized(R)
    w, U = np.linalg.eigh(S)
    accum = 0
    for i in range(n):
        for j in range(n):
            for k in range(n-1):
                prefix = (U[i,k]*U[j,k])**2
                accum += prefix * w[k] * math.exp(2*t*w[k])
    return accum

def get_mutual_information(R, t):
    """
    Get the mutual information between two observations.
    The two observations are of a
    reversible finite-state continuous-time Markov process
    and are separated by time t.
    @param R: rate matrix
    @param t: the amount of time separating the two observations
    """
    return get_expected_ll_ratio(R, t)

def get_mutual_information_diff(R, t):
    # define the number of states
    n = len(R)
    # define the transition matrix and its derivative
    P = scipy.linalg.expm(R*t)
    P_diff = mrate.expm_diff_spectral(R, t)
    # define the stationary distribution
    p = mrate.R_to_distn(R)
    # get the expected log likelihood ratio
    accum = 0
    for i in range(n):
        for j in range(n):
            if p[i] and P[i, j]:
                prefix = p[i] * P_diff[i, j]
                suffix = 1 + math.log(P[i, j]) - math.log(p[j])
                accum += prefix * suffix
    return accum

def get_expected_ll_ratio(R, t):
    """
    This is also the mutual information.
    It is the mutual information between two observations
    of a finite-state continuous-time Markov process at equilibrium
    where the observations are separated by time t.
    """
    # define the number of states
    n = len(R)
    # define the transition matrix
    P = scipy.linalg.expm(R*t)
    # define the stationary distribution
    p = mrate.R_to_distn(R)
    # get the expected log likelihood ratio
    accum = 0
    for i in range(n):
        for j in range(n):
            if p[i] and P[i, j]:
                coeff = p[i] * P[i, j]
                # cancel the p[i] in the numerator and denominator
                #numerator = p[i] * P[i, j]
                #denominator = p[i] * p[j]
                numerator = P[i, j]
                denominator = p[j]
                value = coeff * math.log(numerator / denominator)
                if not np.allclose(np.imag(value), 0):
                    raise ValueError('rogue imaginary number')
                accum += np.real(value)
    return accum


def get_p_id_deriv_ratio(R, t):
    """
    Get (second derivative of p_identity) divided by (first derivative of p_id)
    """
    n = len(R)
    # symmetrize the rate matrix
    v = mrate.R_to_distn(R)
    lam = np.diag(np.sqrt(v))
    rlam = np.diag(np.reciprocal(np.sqrt(v)))
    S = ndot(lam, -R, rlam)
    # eigendecompose the symmetrized rate matrix
    # this should satisfy R = ndot(rlam, V, np.diag(-W), V.T, lam)
    W, V = scipy.linalg.eigh(S)
    # get P and its two derivatives
    P = ndot(rlam, V, np.diag(np.exp(-W*t)), V.T, lam)
    P_dt = ndot(rlam, V, np.diag(-W*np.exp(-W*t)), V.T, lam)
    P_dtt = ndot(rlam, V, np.diag(W*W*np.exp(-W*t)), V.T, lam)
    # get the two derivatives of expected identity
    e_dt = 0.0
    e_dtt = 0.0
    for i in range(n):
        for j in range(n):
            e_dt += v[i] * P_dt[i, i]
            e_dtt += v[i] * P_dtt[i, i]
    return e_dtt / e_dt

def _get_expectation(R, t):
    n = len(R)
    # symmetrize the rate matrix
    v = mrate.R_to_distn(R)
    lam = np.diag(np.sqrt(v))
    rlam = np.diag(np.reciprocal(np.sqrt(v)))
    S = ndot(lam, -R, rlam)
    # eigendecompose the symmetrized rate matrix
    # this should satisfy R = ndot(rlam, V, np.diag(-W), V.T, lam)
    W, V = scipy.linalg.eigh(S)
    # get P and its two derivatives
    P = ndot(rlam, V, np.diag(np.exp(-W*t)), V.T, lam)
    P_dt = ndot(rlam, V, np.diag(-W*np.exp(-W*t)), V.T, lam)
    P_dtt = ndot(rlam, V, np.diag(W*W*np.exp(-W*t)), V.T, lam)
    M = (P*P_dtt - P_dt*P_dt) / P
    expectation = 0.0
    for i in range(n):
        for j in range(n):
            expectation += v[i] * M[i, j]
    return expectation

def _get_expectation_dt(R, t):
    n = len(R)
    # symmetrize the rate matrix
    v = mrate.R_to_distn(R)
    lam = np.diag(np.sqrt(v))
    rlam = np.diag(np.reciprocal(np.sqrt(v)))
    S = ndot(lam, -R, rlam)
    # eigendecompose the symmetrized rate matrix
    # this should satisfy R = ndot(rlam, V, np.diag(-W), V.T, lam)
    W, V = scipy.linalg.eigh(S)
    # get P and its two derivatives
    P = ndot(rlam, V, np.diag(np.exp(-W*t)), V.T, lam)
    P_dt = ndot(rlam, V, np.diag(-W*np.exp(-W*t)), V.T, lam)
    P_dtt = ndot(rlam, V, np.diag(W*W*np.exp(-W*t)), V.T, lam)
    P_dttt = ndot(rlam, V, np.diag(W*W*W*np.exp(-W*t)), V.T, lam)
    M = (P*P*P_dttt + P_dt**3 - 2*P*P_dt*P_dtt) / (P*P)
    expectation_dt = 0.0
    for i in range(n):
        for j in range(n):
            expectation_dt += v[i] * M[i, j]
    return expectation_dt

def get_ml_variance(R, t):
    return -1 / _get_expectation(R, t)

def get_ml_variance_ratio(R, t):
    """
    if v is variance then i get v' / v.
    """
    return -_get_expectation_dt(R, t) / _get_expectation(R, t)


def get_two_state_ml_variance(r, pi0, pi1, t):
    q01 = r*pi1
    q10 = r*pi0
    q00 = -q01
    q11 = -q10
    Q = np.array([
        [q00, q01],
        [q10, q11]])
    P = scipy.linalg.expm(Q*t)
    numerator = P[0,0]*P[1,1]*(P[0,1]+P[1,0])
    denominator = pi0*pi1*r*r*math.exp(-2*r*t)*(1+math.exp(-r*t))
    variance = numerator / denominator
    return variance

class TestDivtime(unittest.TestCase):

    def test_small_variance(self):
        """
        a = .1
        b = .2
        c = .7
        R = np.array([
            [-(b+c), b, c],
            [a, -(a+c), c],
            [a, b, -(a+b)]])
        """
        n = 4
        v = sample_distribution(n)
        S = sample_symmetric_rate_matrix(n)
        R = to_gtr_halpern_bruno(S, v)
        t = 0.0000001
        total_rate = mrate.R_to_total_rate(R)
        var = get_ml_variance(R, t)
        print 'time:', t
        print 'variance:', var
        print 'total rate:', total_rate
        print 'variance per time:', var / t
        print 'reciprocal of total rate:', 1 / total_rate
        print 'total rate times time:', total_rate * t
        print '(reciprocal of total rate) times time:', t / total_rate
        print

    def test_large_variance(self):
        n = 4
        v = sample_distribution(n)
        S = sample_symmetric_rate_matrix(n)
        R = to_gtr_halpern_bruno(S, v)
        """
        a = .1
        b = .2
        c = .7
        R = np.array([
            [-(b+c), b, c],
            [a, -(a+c), c],
            [a, b, -(a+b)]])
        """
        t = 7.5
        dt = 0.0000001
        rtime = mrate.R_to_relaxation_time(R)
        var_a = get_ml_variance(R, t)
        var_b = get_ml_variance(R, t+dt)
        var_slope = (var_b - var_a) / dt
        deriv_ratio = get_p_id_deriv_ratio(R, t)
        clever_ratio = get_ml_variance_ratio(R, t)
        print 'time:', t
        print 'variance:', var_a
        print 'variance slope:', var_slope
        print 'var_slope / var_a:', var_slope / var_a
        print 'var_slope / var_a [clever]:', clever_ratio
        print 'log variance:', math.log(var_a)
        print 'relaxation time:', rtime
        print '2 / relaxation_time:', 2 / rtime
        print "p_id(t)'' / p_id(t)':", deriv_ratio
        print

    def test_variance(self):
        a = .4
        b = .6
        r = 3.0
        R = r * np.array([
            [-a, a],
            [b, -b]])
        t = 0.2
        observed = get_ml_variance(R, t)
        expected = get_two_state_ml_variance(r, a, b, t)
        self.assertTrue(np.allclose(observed, expected))

if __name__ == '__main__':
    unittest.main()

