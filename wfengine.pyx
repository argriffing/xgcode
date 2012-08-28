"""
Wright-Fisher engine.

This is an engine that turns a probability distribution
into a corresponding multinomial distribution
given a population size.
It is vectorized so that it turns multiple probability distributions
into their corresponding multinomial distributions
given a single population size.
The implementation is in Cython for speed
and uses python numpy arrays for speed and convenience.
For compilation instructions see
http://docs.cython.org/src/reference/compilation.html
For example:
$ cython wfengine.pyx
$ gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing \
      -I/usr/include/python2.7 -o wfengine.so wfengine.c
"""

import numpy as np
cimport numpy as np
cimport cython
from libc.math cimport log

np.import_array()

cdef double g_math_neg_inf = float('-inf')

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

@cython.cdivision(True)
cdef double genic_diallelic(int N, int k, double s) nogil:
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

