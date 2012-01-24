"""
Continuous time Markov chain mutual information.

This module is about reversible finite-state
continuous-time Markov processes and the mutual information
between two points in the process separated by a given amount of time.
Reversibility in this context means that the rate matrix
satisfies the detailed balance equations
and implies that the eigenvalues of the rate matrix are real.
This matrix is also assumed to be irreducible.
Combined with reversibility this means that exactly one
of the eigenvalues of the matrix is zero
while all of the other eigenvalues are negative.
Some approximations of the mutual information
have been implemented for very small and very large time separations.
An alternate formulation of mutual information between these random variables
is the expected log likelihood ratio between the joint distribution
of the separated points in the process and the product of their
marginal distributions.
Note that scipy gives eigenvalues in increasing order,
whereas numpy does not make any guarantees about their order.
"""

import math

import numpy as np
import scipy
from scipy import linalg

import mrate

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

def get_mutual_information_small_approx(R, t):
    """
    This is an approximation for small times.
    """
    n = len(R)
    v = mrate.R_to_distn(R)
    S = mrate.symmetrized(R)
    w, U = np.linalg.eigh(S)
    accum = 0
    for i in range(n):
        a = 0
        for k in range(n):
            a += (U[i, k]**2) * math.exp(t * w[k])
        accum += v[i] * a * math.log(a / v[i])
    #print [R[i, i] for i in range(n)]
    #print [sum(U[i, k] * U[i, k] * w[k] for k in range(n)) for i in range(n)]
    #print [sum(U[i, k] * U[i, k] for k in range(n)) for i in range(n)]
    return accum

def get_mutual_information_small_approx_b(R, t):
    """
    This is an approximation for small times.
    Check a decomposition.
    """
    n = len(R)
    v = mrate.R_to_distn(R)
    S = mrate.symmetrized(R)
    w, U = np.linalg.eigh(S)
    accum_a = 0
    accum_b = 0
    accum_c = 0
    accum_d = 0
    for i in range(n):
        a = 0
        b = 0
        for k in range(n):
            prefix = U[i, k] * U[i, k]
            a += prefix * math.exp(t * w[k])
        for k in range(n-1):
            prefix = U[i, k] * U[i, k]
            b += prefix * math.exp(t * w[k])
        x1 = v[i] * v[i]
        x2 = v[i] * b
        y1 = math.log(a)
        y2 = -math.log(v[i])
        accum_a += x1 * y1
        accum_b += x1 * y2
        accum_c += x2 * y1
        accum_d += x2 * y2
    return accum_a + accum_b + accum_c + accum_d

def get_mutual_information_small_approx_c(R, t):
    """
    This is an approximation for small times.
    This is an even more aggressive approximation.
    """
    n = len(R)
    v = mrate.R_to_distn(R)
    S = mrate.symmetrized(R)
    w, U = np.linalg.eigh(S)
    accum = 0
    for i in range(n):
        a = 0
        for k in range(n):
            prefix = U[i, k] * U[i, k]
            a += prefix * math.exp(t * w[k])
        accum += - v[i] * math.log(v[i]) * a
    return accum

def get_mutual_information_small_approx_d(R, t):
    """
    This is an approximation for small times.
    This uses all of the off-diagonal entries of the mutual information
    and also uses an approximation of the off-diagonal entries.
    """
    n = len(R)
    v = mrate.R_to_distn(R)
    S = mrate.symmetrized(R)
    w, U = np.linalg.eigh(S)
    accum_diag_a = 0
    accum_diag_b = 0
    accum_diag_c = 0
    accum_diag_d = 0
    for i in range(n):
        a = 0
        b = 0
        for k in range(n):
            prefix = U[i, k] * U[i, k]
            a += prefix * math.exp(t * w[k])
        for k in range(n-1):
            prefix = U[i, k] * U[i, k]
            b += prefix * math.exp(t * w[k])
        x1 = v[i] * v[i]
        x2 = v[i] * b
        y1 = math.log(a)
        y2 = -math.log(v[i])
        accum_diag_a += x1 * y1
        accum_diag_b += x1 * y2
        accum_diag_c += x2 * y1
        accum_diag_d += x2 * y2
    accum_a = 0
    accum_b = 0
    accum_c = 0
    accum_d = 0
    for i in range(n):
        for j in range(n):
            if i != j:
                prefix = (v[i] * v[j]) ** .5
                a = 0
                for k in range(n):
                    a += U[i, k] * U[j, k] * math.exp(t * w[k])
                b = 0
                for k in range(n-1):
                    b += U[i, k] * U[j, k] * math.exp(t * w[k])
                x1 = v[i] * v[j]
                x2 = prefix * b
                y1 = math.log(a)
                y2 = -math.log(prefix)
                accum_a += x1 * y1
                accum_b += x1 * y2
                accum_c += x2 * y1
                accum_d += x2 * y2
    terms = [
            accum_diag_a, accum_diag_b, accum_diag_c, accum_diag_d,
            accum_a, accum_b, accum_c, accum_d]
    for term in terms:
        print term
    return sum(terms)

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

def get_mutual_information_b(R, t):
    """
    This uses some cancellation.
    """
    n = len(R)
    v = mrate.R_to_distn(R)
    S = mrate.symmetrized(R)
    w, U = np.linalg.eigh(S)
    accum_diag_a = 0
    accum_diag_b = 0
    accum_diag_c = 0
    for i in range(n):
        a = 0
        b = 0
        for k in range(n):
            prefix = U[i, k] * U[i, k]
            a += prefix * math.exp(t * w[k])
        for k in range(n-1):
            prefix = U[i, k] * U[i, k]
            b += prefix * math.exp(t * w[k])
        x1 = v[i] * v[i]
        x2 = v[i] * b
        y1 = math.log(a)
        y2 = -math.log(v[i])
        accum_diag_a += x1 * y1
        accum_diag_b += x1 * y2
        accum_diag_c += x2 * y1
    accum_a = 0
    accum_b = 0
    accum_c = 0
    for i in range(n):
        for j in range(n):
            if i != j:
                prefix = (v[i] * v[j]) ** .5
                a = 0
                for k in range(n):
                    a += U[i, k] * U[j, k] * math.exp(t * w[k])
                b = 0
                for k in range(n-1):
                    b += U[i, k] * U[j, k] * math.exp(t * w[k])
                x1 = v[i] * v[j]
                x2 = prefix * b
                y1 = math.log(a)
                y2 = -math.log(prefix)
                accum_a += x1 * y1
                accum_b += x1 * y2
                accum_c += x2 * y1
    terms = [
            accum_diag_a, accum_diag_b, accum_diag_c,
            accum_a, accum_b, accum_c]
    return sum(terms)

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

def get_mutual_information_diff_b(R, t):
    """
    This is a more symmetrized version.
    Note that two of the three terms are probably structurally zero.
    """
    # get non-spectral summaries
    n = len(R)
    P = scipy.linalg.expm(R*t)
    p = mrate.R_to_distn(R)
    # get spectral summaries
    S = mrate.symmetrized(R)
    w, U = np.linalg.eigh(S)
    G = np.zeros_like(R)
    for i in range(n):
        for j in range(n):
            G[i, j] = 0
            for k in range(n):
                G[i, j] += U[i, k] * U[j, k] * math.exp(t * w[k])
    G_diff = np.zeros_like(R)
    for i in range(n):
        for j in range(n):
            G_diff[i, j] = 0
            for k in range(n):
                G_diff[i, j] += U[i, k] * U[j, k] * w[k] * math.exp(t * w[k])
    B = np.outer(U.T[-1], U.T[-1])
    term_a = np.sum(B * G_diff)
    term_b = np.sum(B * G_diff * np.log(G))
    term_c = -np.sum(B * G_diff * np.log(B))
    #print term_a
    #print term_b
    #print term_c
    return term_b

def get_mutual_information_diff_c(R, t):
    """
    This is a more symmetrized version.
    Some structurally zero terms have been removed.
    """
    # get non-spectral summaries
    n = len(R)
    P = scipy.linalg.expm(R*t)
    p = mrate.R_to_distn(R)
    # get spectral summaries
    S = mrate.symmetrized(R)
    w, U = np.linalg.eigh(S)
    B = np.outer(U.T[-1], U.T[-1])
    G = np.zeros_like(R)
    for i in range(n):
        for j in range(n):
            G[i, j] = 0
            for k in range(n):
                G[i, j] += U[i, k] * U[j, k] * math.exp(t * w[k])
    G_diff = np.zeros_like(R)
    for i in range(n):
        for j in range(n):
            G_diff[i, j] = 0
            for k in range(n):
                G_diff[i, j] += U[i, k] * U[j, k] * w[k] * math.exp(t * w[k])
    return np.sum(B * G_diff * np.log(G))

def get_mutual_information_diff_zero(R):
    """
    Derivative of mutual information at time zero.
    Haha apparently this does not exist.
    """
    # get non-spectral summaries
    n = len(R)
    # get spectral summaries
    S = mrate.symmetrized(R)
    w, U = np.linalg.eigh(S)
    B = np.outer(U.T[-1], U.T[-1])
    G = np.zeros_like(R)
    for i in range(n):
        for j in range(n):
            G[i, j] = 0
            for k in range(n):
                G[i, j] += U[i, k] * U[j, k]
    G_diff = np.zeros_like(R)
    for i in range(n):
        for j in range(n):
            G_diff[i, j] = 0
            for k in range(n):
                G_diff[i, j] += U[i, k] * U[j, k] * w[k]
    print G
    print G_diff
    print B
    return np.sum(B * G_diff * np.log(G))

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

def get_mi_decomposed(U, W, t):
    """
    Get the mutual information at a given time using the decomposition.
    Q = diag(p)^-(1/2) U W U' diag(p)^(1/2)
    The stationary distribution is p.
    Also sqrt(p) is the column of U corresponding to eigenvalue 0.
    """
    pass
