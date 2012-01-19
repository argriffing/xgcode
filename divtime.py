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



############################################################################
# MUTUAL INFORMATION STUFF


def cute_MI_alternate(R, t):
    """
    This is yet another implementation of a large t approximation of MI.
    It is related the expectation of the ratio of the probability
    of what you actually saw to the probability of seeing
    what you saw given independence.
    It is half of one less than this expectation.
    It is not as numerically stable as other large t approximations.
    """
    # define the number of states
    n = len(R)
    # define the transition matrix
    P = scipy.linalg.expm(R*t)
    # define the stationary distribution
    p = mrate.R_to_distn(R)
    s = np.sqrt(p)
    # get the expected log likelihood ratio
    accum = 0
    for i in range(n):
        for j in range(n):
            p_joint = p[i] * P[i, j]
            p_independent = p[i] * p[j]
            accum += p_joint * (p_joint / p_independent)
    return (accum - 1) / 2

def cute_MI_alternate_b(R, t):
    """
    It should closely approximate mutual information when t is not tiny.
    """
    # define the number of states
    n = len(R)
    # define the transition matrix
    P = scipy.linalg.expm(R*t)
    # define the stationary distribution
    p = mrate.R_to_distn(R)
    s = np.sqrt(p)
    # get the expected log likelihood ratio
    accum = 0
    for i in range(n):
        for j in range(n):
            p_joint = p[i] * P[i, j]
            value = p_joint / (s[i] * s[j]) - (s[i] * s[j])
            accum += (value * value) / 2
    return accum

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

def get_mutual_information_approx_c(R, t):
    """
    This is an approximation for large times.
    It has been rewritten using orthogonality.
    It has also been rewritten using orthonormality.
    """
    n = len(R)
    v = mrate.R_to_distn(R)
    S = mrate.symmetrized(R)
    w, U = np.linalg.eigh(S)
    accum = 0
    for k in range(n-1):
        accum += math.exp(2*t*w[k])
    return accum / 2

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

def get_mutual_information_diff_approx_c(R, t):
    """
    This is an approximation for large times.
    It has been rewritten using orthogonality.
    It has also been rewritten using orthonormality.
    """
    n = len(R)
    v = mrate.R_to_distn(R)
    S = mrate.symmetrized(R)
    w, U = np.linalg.eigh(S)
    accum = 0
    for k in range(n-1):
        accum += w[k]*math.exp(2*t*w[k])
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


############################################################################
# ASYMPTOTIC VARIANCE STUFF

def get_asymptotic_variance(R, t):
    """
    Asymptotic variance is the negative reciprocal of an expectation.
    The expectation is of the second derivative of the log likelihood.
    Use a fact about the second derivative of logarithm
    to evaluate this asymptotic variance.
    Also this is probably a dumb non-spectral way that will be improved later.
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
            # define the second derivative of f
            f_dtt = 0
            for k in range(n):
                f_dtt += U[i, k] * U[j, k] * w[k] * w[k] * math.exp(t * w[k])
            f_dtt *= (p[i] * p[j])**.5
            # accumulate the contribution of this entry to the expectation
            accum += f_dtt - (f_dt * f_dt) / f
    return - 1 / accum

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
        print 'mutual information:', get_mutual_information(R, t)
        print 'reciprocal of MI:', 1.0 / get_mutual_information(R, t)
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

