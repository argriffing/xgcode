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
    This is deliberately not a python function,
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
        np.ndarray[np.float64_t, ndim=1] log_factorial,
        ):
    """
    @param N: population size
    @param k: number of bins
    @param log_distn: the log of distribution over bins
    @param counts: the observed counts over bins
    @param log_factorial: precomputed logarithm of factorial values
    @return: logarithm of pmf
    """
    cdef int i
    cdef int c
    cdef double accum = log_factorial[N]
    # check for an impossible state
    for i in range(k):
        if counts[i] and log_distn[i] == g_math_neg_inf:
            return g_math_neg_inf
    # add the contribution of the counts to the multinomial coefficient
    # add the contribution of probabilities
    for i in range(k):
        c = counts[i]
        accum -= log_factorial[c]
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
    log_factorial = np.zeros(N+1)
    accum = 0
    for i in range(2, N+1):
        accum += log(i)
        log_factorial[i] = accum
    # Make the array.
    for i, compo_list in enumerate(gen_population_compositions(N, k)):
        compo = np.array(compo_list, dtype=np.int)
        for j in range(ndistns):
            # Compute the log multinomial probability.
            M[j, i] = multinomial_log_pmf(
                    N, k, log_distns[j], compo, log_factorial)
    # Return the array.
    return M

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
    cdef np.ndarray[np.float64_t, ndim=1] log_factorial = np.zeros(N+1)
    accum = 0
    for i in range(2, N+1):
        accum += log(i)
        log_factorial[i] = accum
    # Make the array.
    cdef np.ndarray[np.int_t, ndim=2] compos
    compos = np.array(list(gen_population_compositions(N, k)))
    for i in range(ncols_out):
        # get the composition as a numpy array
        #compo = np.array(compo_list, dtype=np.int)
        # define the log of multinomial coefficient for this composition
        log_multinomial_coeff = log_factorial[N]
        for index in range(k):
            log_multinomial_coeff -= log_factorial[compos[i, index]]
        for j in range(ndistns):
            # Compute the log multinomial probability.
            accum = log_multinomial_coeff
            for index in range(k):
                if compos[i, index]:
                    accum += compos[i, index] * log_distns[j, index]
            M[j, i] = accum
    # Return the array.
    return M

