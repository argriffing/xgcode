import unittest
import math

import numpy as np
import scipy
from scipy import linalg

import mrate
import ctmcmi
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


############################################################################
# ASYMPTOTIC VARIANCE STUFF


def get_fisher_information(R, t):
    # get non-spectral summaries
    n = len(R)
    P = scipy.linalg.expm(R*t)
    p = mrate.R_to_distn(R)
    # get spectral summaries
    S = mrate.symmetrized(R)
    w, U = scipy.linalg.eigh(S)
    # compute the asymptotic variance
    accum = 0
    for i in range(n):
        for j in range(n):
            # define f
            f = p[i] * P[i, j]
            # define the first derivative of f
            f_dt = 0
            for k in range(n):
                f_dt += U[i, k] * U[j, k] * w[k] * math.exp(t * w[k])
            f_dt *= (p[i] * p[j])**.5
            # define the second derivative of f
            f_dtt = 0
            for k in range(n):
                f_dtt += U[i, k] * U[j, k] * w[k] * w[k] * math.exp(t * w[k])
            f_dtt *= (p[i] * p[j])**.5
            # accumulate the contribution of this entry to the expectation
            accum += f_dtt - (f_dt * f_dt) / f
    return -accum


def get_asymptotic_variance(R, t):
    """
    Asymptotic variance is the negative reciprocal of an expectation.
    The expectation is of the second derivative of the log likelihood.
    Use a fact about the second derivative of logarithm
    to evaluate this asymptotic variance.
    Also this is probably a dumb non-spectral way that will be improved later.
    """
    return 1 / get_fisher_information(R, t)

def get_asymptotic_variance_b(R, t):
    """
    Break up the sum into two parts and investigate each separately.
    The second part with only the second derivative is zero.
    """
    # get non-spectral summaries
    n = len(R)
    P = scipy.linalg.expm(R*t)
    p = mrate.R_to_distn(R)
    # get spectral summaries
    S = mrate.symmetrized(R)
    w, U = np.linalg.eigh(S)
    # compute the asymptotic variance
    accum_a = 0
    for i in range(n):
        for j in range(n):
            # define f
            f = p[i] * P[i, j]
            # define the first derivative of f
            f_dt = 0
            for k in range(n):
                f_dt += U[i, k] * U[j, k] * w[k] * math.exp(t * w[k])
            f_dt *= (p[i] * p[j])**.5
            accum_a -= (f_dt * f_dt) / f
    accum_b = 0
    for i in range(n):
        for j in range(n):
            # define the second derivative of f
            f_dtt = 0
            for k in range(n):
                f_dtt += U[i, k] * U[j, k] * w[k] * w[k] * math.exp(t * w[k])
            f_dtt *= (p[i] * p[j])**.5
            # accumulate the contribution of this entry to the expectation
            accum_b += f_dtt
    return - 1 / (accum_a + accum_b)

def get_asymptotic_variance_c(R, t):
    """
    Re-evaluate, this time throwing away the part that is structurally zero.
    """
    # get non-spectral summaries
    n = len(R)
    P = scipy.linalg.expm(R*t)
    p = mrate.R_to_distn(R)
    # get spectral summaries
    S = mrate.symmetrized(R)
    w, U = np.linalg.eigh(S)
    # compute the asymptotic variance
    accum = 0
    for i in range(n):
        for j in range(n):
            # define f
            f = p[i] * P[i, j]
            # define the first derivative of f
            f_dt = 0
            for k in range(n):
                f_dt += U[i, k] * U[j, k] * w[k] * math.exp(t * w[k])
            f_dt *= (p[i] * p[j])**.5
            accum -= (f_dt * f_dt) / f
    return - 1 / accum

def get_asymptotic_variance_d(R, t):
    """
    Use a very aggressive approximation.
    """
    # get non-spectral summaries
    n = len(R)
    P = scipy.linalg.expm(R*t)
    p = mrate.R_to_distn(R)
    # get spectral summaries
    S = mrate.symmetrized(R)
    w, U = np.linalg.eigh(S)
    # compute the asymptotic variance approximation
    accum = 0
    for k in range(n-1):
        accum += w[k] * w[k] * math.exp(2 * t * w[k])
    return 1 / accum

def get_asymptotic_variance_e(R, t):
    """
    Try to mitigate the damage of the aggressive approximation.
    The next step is to try to simplify this complicated correction.
    But I have not been able to do this.
    """
    # get non-spectral summaries
    n = len(R)
    P = scipy.linalg.expm(R*t)
    p = mrate.R_to_distn(R)
    # get spectral summaries
    S = mrate.symmetrized(R)
    w, U = np.linalg.eigh(S)
    # compute the asymptotic variance approximation
    accum = 0
    for k in range(n-1):
        accum += w[k] * w[k] * math.exp(2 * t * w[k])
    accum_b = 0
    G_a = np.zeros_like(R)
    G_b = np.zeros_like(R)
    for i in range(n):
        for j in range(n):
            prefix = (p[i] * p[j]) ** -.5
            a = 0
            for k in range(n-1):
                a += U[i, k] * U[j, k] * math.exp(t * w[k])
            b = 0
            for k in range(n-1):
                b += U[i, k] * U[j, k] * w[k] * math.exp(t * w[k])
            suffix = a * b * b
            value = prefix * suffix
            accum_b += value
    return 1 / (accum - accum_b)


############################################################################
# UNKNOWN OLD STUFF



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
    return - 1 / _get_expectation(R, t)

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
        R = mrate.to_gtr_halpern_bruno(S, v)
        t = 0.0000001
        total_rate = mrate.Q_to_expected_rate(R)
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
        R = mrate.to_gtr_halpern_bruno(S, v)
        """
        a = .1
        b = .2
        c = .7
        R = np.array([
            [-(b+c), b, c],
            [a, -(a+c), c],
            [a, b, -(a+b)]])
        """
        t = 5.0
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
        print '--- new attempt ---'
        print 'mutual information:', ctmcmi.get_mutual_information(R, t)
        print 'reciprocal of MI:', 1.0 / ctmcmi.get_mutual_information(R, t)
        print 'asymptotic variance:', get_asymptotic_variance(R, t)
        print 'asymptotic variance (ver. 2):', get_asymptotic_variance_b(R, t)
        print 'asymptotic variance (ver. 3):', get_asymptotic_variance_c(R, t)
        print 'AV approx (ver. 4):', get_asymptotic_variance_d(R, t)
        print 'AV approx (ver. 5):', get_asymptotic_variance_e(R, t)
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
    np.set_printoptions(linewidth=200)
    unittest.main()

