"""
Wright-Fisher engine.

This is an engine that turns a probability distribution
into a corresponding multinomial distribution
given a population size.
It is vectorized so that it turns multiple probability distributions
into their corresponding multinomial distributions
given a single population size.
.
The implementation is in Cython for speed
and uses python numpy arrays for speed and convenience.
For compilation instructions see
http://docs.cython.org/src/reference/compilation.html
For example:
$ cython wfengine.pyx
$ gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing \
      -I/usr/include/python2.7 -o wfengine.so wfengine.c
"""

# TODO generalize the diallelic transition matrix creation functions

import math
import numpy as np
cimport numpy as np
cimport cython
from libc.math cimport log, exp

np.import_array()

cdef double g_math_neg_inf = float('-inf')

def binomial_coefficient(n, k):
    """
    Modified from a function by Andrew Dalke.
    This is deliberately a python function rather than a cdef,
    because we want to allow large integers.
    """
    if 0 <= k <= n:
        ntok = 1
        ktok = 1
        for t in range(1, min(k, n - k) + 1):
            ntok *= n
            ktok *= t
            n -= 1
        return ntok // ktok
    else:
        return 0

#XXX useful but unused
def invert_binomial_coefficient(n_choose_k, k):
    """
    This is a plain python function and works with large integers.
    It extracts n given n_choose_k and k.
    The complexity is linear with k.
    """
    if k < 1:
        raise ValueError('k should be at least 1')
    n_low = (math.factorial(k) * n_choose_k)**(1.0 / k)
    n_high = n_low + k
    for n in range(int(math.floor(n_low)), int(math.ceil(n_high)) + 1):
        if binomial_coefficient(n, k) == n_choose_k:
            return n
    raise ValueError('failed to invert binomial coefficient')

@cython.boundscheck(False)
@cython.wraparound(False)
def get_lmcs(
        np.ndarray[np.int_t, ndim=2] M,
        ):
    """
    Logs of multinomial coefficients.
    @param M: M[i,j] is count of bin j in state i
    @return: a one dimensional array of log multinomial coefficients
    """
    cdef int nstates = M.shape[0]
    cdef int k = M.shape[1]
    cdef int N = np.sum(M[0])
    cdef int i
    cdef np.ndarray[np.float64_t, ndim=1] log_fact = get_log_fact_array(N)
    cdef np.ndarray[np.float64_t, ndim=1] v = np.empty(nstates)
    for i in range(nstates):
        v[i] = log_fact[N]
        for index in range(k):
            v[i] -= log_fact[M[i, index]]
    return v

#XXX deprecated
def gen_population_compositions(int N, int k):
    """
    Yield (N+k-1 choose N) compositions of length k.
    @param N: population size
    @param k: number of bins
    """
    if k == 1:
        yield [N]
    elif k == 2:
        for i in range(N+1):
            yield [i, N-i]
    else:
        for i in range(N+1):
            for suffix in gen_population_compositions(N-i, k-1):
                yield [i] + suffix

# This is from a cython tutorial.
cdef inline int int_min(int a, int b):
    return a if a >= b else b

@cython.boundscheck(False)
@cython.wraparound(False)
cdef double multinomial_log_pmf(
        int N,
        int k,
        np.ndarray[np.float64_t, ndim=1] log_distn,
        np.ndarray[np.int_t, ndim=1] counts,
        np.ndarray[np.float64_t, ndim=1] log_fact,
        ):
    """
    @param N: population size
    @param k: number of bins
    @param log_distn: the log of distribution over bins
    @param counts: the observed counts over bins
    @param log_fact: precomputed logarithm of factorial values
    @return: logarithm of pmf
    """
    cdef int i
    cdef int c
    cdef double accum = log_fact[N]
    # check for an impossible state
    for i in range(k):
        if counts[i] and log_distn[i] == g_math_neg_inf:
            return g_math_neg_inf
    # add the contribution of the counts to the multinomial coefficient
    # add the contribution of probabilities
    for i in range(k):
        c = counts[i]
        accum -= log_fact[c]
        if c:
            accum += c * log_distn[i]
    return accum

def expand_multinomials_slow(N, log_distns):
    """
    Each row of distns is a distribution over k bins.
    Each row of the returned array
    is a multinomial distribution over the ways to put N things into k bins.
    The returned array has the same number of columns as the distn array,
    and the entries are logarithms of probabilities.
    @param N: integer population size
    @param log_distns: numpy 2d array with ndistns rows and k columns
    @return: numpy 2d array with ndistns rows and choose(N+k-1,N) columns
    """
    ndistns, k = log_distns.shape
    # Allocate the output numpy array
    # whose rows will have the multinomial distributions.
    ncols_out = binomial_coefficient(N+k-1, N)
    M = np.zeros((ndistns, ncols_out))
    # Precompute some logarithms of factorials up to N.
    log_fact = np.zeros(N+1)
    accum = 0
    for i in range(2, N+1):
        accum += log(i)
        log_fact[i] = accum
    # Make the array.
    for i, compo_list in enumerate(gen_population_compositions(N, k)):
        compo = np.array(compo_list, dtype=np.int)
        for j in range(ndistns):
            # Compute the log multinomial probability.
            M[j, i] = multinomial_log_pmf(
                    N, k, log_distns[j], compo, log_fact)
    # Return the array.
    return M

@cython.boundscheck(False)
@cython.wraparound(False)
cdef get_log_fact_array(int N):
    """
    Precompute some logarithms of factorials up to N.
    @param N: max integer whose log of factorial to compute
    @return: a numpy array of length N+1
    """
    cdef np.ndarray[np.float64_t, ndim=1] log_fact = np.zeros(N+1)
    cdef double accum = 0
    for i in range(2, N+1):
        accum += log(i)
        log_fact[i] = accum
    return log_fact

@cython.boundscheck(False)
@cython.wraparound(False)
def expand_multinomials(
        int N,
        np.ndarray[np.float64_t, ndim=2] log_distns,
        ):
    """
    Each row of distns is a distribution over k bins.
    Each row of the returned array
    is a multinomial distribution over the ways to put N things into k bins.
    The returned array has the same number of columns as the distn array,
    and the entries are logarithms of probabilities.
    @param N: integer population size
    @param log_distns: numpy 2d array with ndistns rows and k columns
    @return: numpy 2d array with ndistns rows and choose(N+k-1,N) columns
    """
    cdef int i, j, k, index
    cdef int ndistns, ncolsout
    cdef double accum, log_multinomial_coeff
    # initialize
    ndistns = log_distns.shape[0]
    k = log_distns.shape[1]
    # Allocate the output numpy array
    # whose rows will have the multinomial distributions.
    ncols_out = binomial_coefficient(N+k-1, N)
    cdef np.ndarray[np.float64_t, ndim=2] M = np.zeros((ndistns, ncols_out))
    # Precompute some logarithms of factorials up to N.
    cdef np.ndarray[np.float64_t, ndim=1] log_fact = get_log_fact_array(N)
    # Make the array.
    cdef np.ndarray[np.int_t, ndim=2] compos
    compos = np.array(list(gen_population_compositions(N, k)))
    for i in range(ncols_out):
        # define the log of multinomial coefficient for this composition
        log_multinomial_coeff = log_fact[N]
        for index in range(k):
            log_multinomial_coeff -= log_fact[compos[i, index]]
        for j in range(ndistns):
            # Compute the log multinomial probability.
            accum = log_multinomial_coeff
            for index in range(k):
                if compos[i, index]:
                    accum += compos[i, index] * log_distns[j, index]
            M[j, i] = accum
    # Return the array.
    return M

@cython.boundscheck(False)
@cython.wraparound(False)
cdef mvhyperg(
        np.ndarray[np.int_t, ndim=1] s_large,
        np.ndarray[np.int_t, ndim=1] s_small,
        np.ndarray[np.float64_t, ndim=1] log_fact,
        ):
    """
    Multivariate hypergeometric log likelihood.
    This is named similarly to the fortran function in flib.
    @return: log likelihood
    """
    # get the number of bins
    cdef int k = s_large.shape[0]
    cdef int i
    # check for impossible selections
    for i in range(k):
        if s_small[i] > s_large[i]:
            return g_math_neg_inf
    # get the number of balls
    cdef int N = 0
    cdef int n = 0
    cdef int v_N, v_n
    cdef double log_num = 0
    cdef double log_den = 0
    for i in range(k):
        v_N = s_large[i]
        v_n = s_small[i]
        N += v_N
        n += v_n
        log_num += log_fact[v_N] - log_fact[v_n] - log_fact[v_N - v_n]
    log_den += log_fact[N] - log_fact[n] - log_fact[N-n]
    return log_num - log_den

@cython.boundscheck(False)
@cython.wraparound(False)
def reduce_hypergeometric(
        np.ndarray[np.int_t, ndim=2] M_large,
        np.ndarray[np.int_t, ndim=2] M_small,
        np.ndarray[np.float64_t, ndim=1] w_large,
        ):
    """
    @param M_large: states defined by N balls in k bins
    @param M_small: states defined by n balls in k bins
    @param w_large: weight per large state
    @return: weight per small state
    """
    # nstates_large is the number of larger samples
    # nstates_small is the number of smaller samples
    # k is the number of unique microstates.
    # N is the size of the larger sample.
    # n is the size of the smaller sample.
    cdef int nstates_large = M_large.shape[0]
    cdef int nstates_small = M_small.shape[0]
    cdef int k = M_large.shape[1]
    cdef int N = M_large[0].sum()
    cdef int n = M_small[0].sum()
    cdef int i, j
    cdef double loglik
    cdef np.ndarray[np.float64_t, ndim=1] w_small = np.zeros(nstates_small)
    cdef np.ndarray[np.float64_t, ndim=1] log_fact = get_log_fact_array(N)
    for i in range(nstates_small):
        for j in range(nstates_large):
            loglik = mvhyperg(M_large[j], M_small[i], log_fact)
            w_small[i] += w_large[j] * exp(loglik)
    return w_small

@cython.cdivision(True)
cdef double diallelic_chen(
        int N, int k, double fAA, double faA, double faa) nogil:
    """
    Use the notation of Christina Chen et al.
    The title of the paper is
    Effects of dominance on the probability of fixation of a mutant allele.
    """
    cdef int r = N - k
    cdef double AA = fAA*k*k
    cdef double aA = faA*k*r
    cdef double aa = faa*r*r
    return (AA + aA) / (AA + 2*aA + aa)

@cython.cdivision(True)
cdef double diallelic_recessive(int N, int k, double s, double h) nogil:
    """
    This uses selection notation analogous to Kai Zeng 2010.
    When 0 < s <= 1 the first allele is more fit.
    When 0 <= h < 1/2 the first allele is recessive.
    When h = 1/2 the selection is genic additive.
    When 1/2 < h <= 1 the first allele is dominant.
    """
    cdef int r = N - k
    cdef double f00 = 1.0
    cdef double f11 = 1.0 - s
    cdef double f01 = h * f00 + (1-h) * f11
    cdef double aa = f00*k*k
    cdef double ab = f01*k*r
    cdef double bb = f11*r*r
    return (aa + ab) / (aa + 2*ab + bb)

@cython.cdivision(True)
cdef double genic_diallelic(int N, int k, double s) nogil:
    """
    This uses the notation in Kai Zeng 2010.
    """
    cdef int r = N - k
    cdef double f00 = 1.0
    cdef double f11 = 1.0 - s
    cdef double f01 = 0.5 * (f00 + f11)
    cdef double aa = f00*k*k
    cdef double ab = f01*k*r
    cdef double bb = f11*r*r
    return (aa + ab) / (aa + 2*ab + bb)

@cython.cdivision(True)
cdef double genic_diallelic_ohta(int N, int k, double s) nogil:
    """
    This corresponds to using 1+s and 1 as opposed to using 1 and 1-s.
    """
    cdef double p = k / (1.0 * N)
    cdef double delta = (0.5 * s) * p * (1 - p) / (1 + s*p)
    return p + delta

@cython.boundscheck(False)
@cython.wraparound(False)
def create_genic_diallelic(int N_diploid, double s):
    """
    Create a Wright-Fisher transition matrix.
    Use genic selection with two alleles and a diploid population.
    The returned transition matrix will be a square matrix as a numpy array,
    and it will have 2N+1 states.
    State k corresponds to the presence of k preferred alleles
    in the population when s is positive.
    @param N_diploid: diploid population size
    @param s: an additive selection value that is positive by convention
    @return: entrywise logarithms of a transition matrix
    """
    cdef int N = N_diploid * 2
    # declare intermediate variables
    cdef int i, j
    cdef double p, log_p, log_pcompl
    # init the transition matrix
    cdef np.ndarray[np.float64_t, ndim=2] M = np.zeros((N+1,  N+1))
    # Precompute some logarithms of factorials up to N.
    cdef np.ndarray[np.float64_t, ndim=1] log_fact = get_log_fact_array(N)
    # i and j are the number of preferred alleles
    # in the parent and child generations respectively
    for i in range(N+1):
        p = genic_diallelic(N, i, s)
        log_p = log(p)
        log_pcompl = log(1-p)
        for j in range(N+1):
            M[i, j] = log_fact[N] - log_fact[j] - log_fact[N-j]
            if j:
                M[i, j] += j * log_p
            if N-j:
                M[i, j] += (N-j) * log_pcompl
    return M

@cython.boundscheck(False)
@cython.wraparound(False)
def create_genic_diallelic_ohta(int N_diploid, double s):
    # TODO this is mostly copy and paste
    cdef int N = N_diploid * 2
    # declare intermediate variables
    cdef int i, j
    cdef double p, log_p, log_pcompl
    # init the transition matrix
    cdef np.ndarray[np.float64_t, ndim=2] M = np.zeros((N+1,  N+1))
    # Precompute some logarithms of factorials up to N.
    cdef np.ndarray[np.float64_t, ndim=1] log_fact = get_log_fact_array(N)
    # i and j are the number of preferred alleles
    # in the parent and child generations respectively
    for i in range(N+1):
        p = genic_diallelic_ohta(N, i, s)
        log_p = log(p)
        log_pcompl = log(1-p)
        for j in range(N+1):
            M[i, j] = log_fact[N] - log_fact[j] - log_fact[N-j]
            if j:
                M[i, j] += j * log_p
            if N-j:
                M[i, j] += (N-j) * log_pcompl
    return M

@cython.boundscheck(False)
@cython.wraparound(False)
def create_diallelic_recessive(int N_diploid, double s, double h):
    # TODO this is mostly copy and paste
    cdef int N = N_diploid * 2
    # declare intermediate variables
    cdef int i, j
    cdef double p, log_p, log_pcompl
    # init the transition matrix
    cdef np.ndarray[np.float64_t, ndim=2] M = np.zeros((N+1,  N+1))
    # Precompute some logarithms of factorials up to N.
    cdef np.ndarray[np.float64_t, ndim=1] log_fact = get_log_fact_array(N)
    # i and j are the number of preferred alleles
    # in the parent and child generations respectively
    for i in range(N+1):
        p = diallelic_recessive(N, i, s, h)
        log_p = log(p)
        log_pcompl = log(1-p)
        for j in range(N+1):
            M[i, j] = log_fact[N] - log_fact[j] - log_fact[N-j]
            if j:
                M[i, j] += j * log_p
            if N-j:
                M[i, j] += (N-j) * log_pcompl
    return M

@cython.boundscheck(False)
@cython.wraparound(False)
def create_diallelic_chen(int N_diploid, double fAA, double faA, double faa):
    # TODO this is mostly copy and paste
    cdef int N = N_diploid * 2
    # declare intermediate variables
    cdef int i, j
    cdef double p, log_p, log_pcompl
    # init the transition matrix
    cdef np.ndarray[np.float64_t, ndim=2] M = np.zeros((N+1,  N+1))
    # Precompute some logarithms of factorials up to N.
    cdef np.ndarray[np.float64_t, ndim=1] log_fact = get_log_fact_array(N)
    # i and j are the number of preferred alleles
    # in the parent and child generations respectively
    for i in range(N+1):
        p = diallelic_chen(N, i, fAA, faA, faa)
        log_p = log(p)
        log_pcompl = log(1-p)
        for j in range(N+1):
            M[i, j] = log_fact[N] - log_fact[j] - log_fact[N-j]
            if j:
                M[i, j] += j * log_p
            if N-j:
                M[i, j] += (N-j) * log_pcompl
    return M

@cython.boundscheck(False)
@cython.wraparound(False)
def create_genic(
        np.ndarray[np.float64_t, ndim=1] lmcs,
        np.ndarray[np.float64_t, ndim=2] lps,
        np.ndarray[np.int_t, ndim=2] M,
        ):
    """
    This is a more flexible way to create a multinomial transition matrix.
    It is also fast.
    The output may have -inf but it should not have nan.
    @param lmcs: log multinomial count per state
    @param lps: log probability per haplotype per state
    @param M: allele count per haplotype per state
    @return: entrywise log of transition matrix
    """
    cdef int nstates = M.shape[0]
    cdef int k = M.shape[1]
    cdef int i, j, index
    cdef np.ndarray[np.float64_t, ndim=2] L = np.zeros((nstates, nstates))
    for i in range(nstates):
        for j in range(nstates):
            L[i, j] = lmcs[j]
            for index in range(k):
                if M[j, index]:
                    L[i, j] += M[j, index] * lps[i, index]
    return L

