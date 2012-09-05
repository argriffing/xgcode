"""
Quickly create large mutation and recombination rate matrices.

This is for studying compensatory nucleotide changes in evolution.
The four haplotypes are (AB, Ab, aB, ab).
The AB and ab haplotypes have high fitness,
while the Ab and aB haplotypes have low fitness.
.
The implementation is in Cython for speed
and uses python numpy arrays for speed and convenience.
For compilation instructions see
http://docs.cython.org/src/reference/compilation.html
For example:
$ cython wfcompens.pyx
$ gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing \
      -I/usr/include/python2.7 -o wfcompens.so wfcompens.c
"""

import numpy as np
cimport numpy as np
cimport cython
from libc.math cimport log

np.import_array()

cdef int AB_type = 0
cdef int Ab_type = 1
cdef int aB_type = 2
cdef int ab_type = 3

@cython.boundscheck(False)
@cython.wraparound(False)
def create_mutation(
        np.ndarray[np.int_t, ndim=2] M,
        np.ndarray[np.int_t, ndim=4] T,
        ):
    """
    The scaling of the resulting rate matrix is strange.
    Every rate is an integer, but in double precision float format.
    @param M: M[i,j] is the count of allele j in state index i
    @param T: T[AA, Ab, aB, ab] is the state index for the given allele counts
    @return: mutation rate matrix
    """
    cdef int i
    cdef int AB, Ab, aB, ab
    cdef int nstates = M.shape[0]
    cdef np.ndarray[np.float64_t, ndim=2] R = np.zeros((nstates, nstates))
    for i in range(nstates):
        # meticulously unpack the allele counts from the corresponding state
        AB = M[i, 0]
        Ab = M[i, 1]
        aB = M[i, 2]
        ab = M[i, 3]
        #
        if AB > 0:
            R[i, T[AB-1, Ab+1, aB,   ab  ]] = AB
            R[i, T[AB-1, Ab,   aB+1, ab  ]] = AB
        if Ab > 0:
            R[i, T[AB+1, Ab-1, aB,   ab  ]] = Ab
            R[i, T[AB,   Ab-1, aB,   ab+1]] = Ab
        if aB > 0:
            R[i, T[AB+1, Ab,   aB-1, ab  ]] = aB
            R[i, T[AB,   Ab,   aB-1, ab+1]] = aB
        if ab > 0:
            R[i, T[AB,   Ab+1, aB,   ab-1]] = ab
            R[i, T[AB,   Ab,   aB+1, ab-1]] = ab
        R[i, i] = -2*(AB + Ab + aB + ab)
    return R

@cython.boundscheck(False)
@cython.wraparound(False)
def create_recomb(
        np.ndarray[np.int_t, ndim=2] M,
        np.ndarray[np.int_t, ndim=4] T,
        ):
    """
    The scaling of the resulting rate matrix is strange.
    Every rate is an integer, but in double precision float format.
    @param M: M[i,j] is the count of allele j in state index i
    @param T: T[AA, Ab, aB, ab] is the state index for the given allele counts
    @return: recombination rate matrix
    """
    cdef int i
    cdef int AB, Ab, aB, ab
    cdef int AB_ab, Ab_aB
    cdef int nstates = M.shape[0]
    cdef np.ndarray[np.float64_t, ndim=2] R = np.zeros((nstates, nstates))
    for i in range(nstates):
        # meticulously unpack the allele counts from the corresponding state
        AB = M[i, 0]
        Ab = M[i, 1]
        aB = M[i, 2]
        ab = M[i, 3]
        #
        AB_ab = AB * ab
        Ab_aB = Ab * aB
        #
        if AB_ab > 0:
            R[i, T[AB-1, Ab+1, aB+1, ab-1]] = AB_ab
        if Ab_aB > 0:
            R[i, T[AB+1, Ab-1, aB-1, ab+1]] = Ab_aB
        R[i, i] = -(AB_ab + Ab_aB)
    return R

@cython.boundscheck(False)
@cython.wraparound(False)
def create_selection(
        double s,
        np.ndarray[np.int_t, ndim=2] M,
        ):
    if s >= 1:
        raise ValueError(
                'selection s must be less than 1 '
                'but observed: %s' % s)
    #
    cdef double AB, Ab, aB, ab
    #
    cdef int nstates = M.shape[0]
    cdef int k = M.shape[1]
    cdef double neg_logp
    cdef np.ndarray[np.float64_t, ndim=2] L = np.empty((nstates, k))
    for i in range(nstates):
        AB = M[i, 0]
        Ab = M[i, 1]
        aB = M[i, 2]
        ab = M[i, 3]
        neg_logp = -log(AB + ab + (1-s)*(Ab + aB))
        L[i, 0] = neg_logp + log(AB)
        L[i, 1] = neg_logp + log(Ab*(1-s))
        L[i, 2] = neg_logp + log(aB*(1-s))
        L[i, 3] = neg_logp + log(ab)
    return L



###########################################################################
# Forward simulation helper functions.
# These are intended for speed but may be premature optimization.
# Faster speeds could be achieved through gsl random sample generation.


@cython.boundscheck(False)
@cython.wraparound(False)
def multiple_mutation(
        np.ndarray[np.int_t, ndim=1] mutable_state,
        np.ndarray[np.int_t, ndim=1] mutable_counts,
        np.ndarray[np.int_t, ndim=1] sampled_individuals,
        np.ndarray[np.int_t, ndim=1] sampled_loci,
        ):
    """
    All of the random sampling has already occurred.
    """
    cdef int nsamples = sampled_individuals.shape[0]
    cdef int i, index
    cdef int old_type, new_type
    cdef int locus
    for i in range(nsamples):
        index = sampled_individuals[i]
        locus = sampled_loci[i]
        old_type = mutable_state[index]
        if old_type == AB_type:
            if not locus:
                new_type = aB_type
            else:
                new_type = Ab_type
        elif old_type == Ab_type:
            if not locus:
                new_type = ab_type
            else:
                new_type = AB_type
        elif old_type == aB_type:
            if not locus:
                new_type = AB_type
            else:
                new_type = ab_type
        elif old_type == ab_type:
            if not locus:
                new_type = Ab_type
            else:
                new_type = aB_type
        mutable_state[index] = new_type
        mutable_counts[old_type] -= 1
        mutable_counts[new_type] += 1
    return None

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def multiple_recombination(
        np.ndarray[np.int_t, ndim=1] mutable_state,
        np.ndarray[np.int_t, ndim=1] mutable_counts,
        np.ndarray[np.int_t, ndim=1] sampled_individuals,
        ):
    """
    All of the random sampling has already occurred.
    @param sampled_individuals: should be of even length
    """
    cdef int nsamples = sampled_individuals.shape[0] / 2
    cdef int i
    cdef int index_0, index_1
    cdef int t0, t1
    for i in range(nsamples):
        index_0 = sampled_individuals[i*2]
        index_1 = sampled_individuals[i*2 + 1]
        t0 = mutable_state[index_0]
        t1 = mutable_state[index_1]
        if (t0 == AB_type and t1 == ab_type) or (
                t0 == ab_type and t1 == AB_type):
            mutable_state[index_0] = aB_type
            mutable_state[index_1] = Ab_type
            mutable_counts[AB_type] -= 1
            mutable_counts[ab_type] -= 1
            mutable_counts[Ab_type] += 1
            mutable_counts[aB_type] += 1
        elif (t0 == Ab_type and t1 == aB_type) or (
                t0 == aB_type and t1 == Ab_type):
            mutable_state[index_0] = AB_type
            mutable_state[index_1] = ab_type
            mutable_counts[AB_type] += 1
            mutable_counts[ab_type] += 1
            mutable_counts[Ab_type] -= 1
            mutable_counts[aB_type] -= 1
    return None

@cython.boundscheck(False)
@cython.wraparound(False)
def expand_counts(
        np.ndarray[np.int_t, ndim=1] mutable_state,
        np.ndarray[np.int_t, ndim=1] counts,
        ):
    cdef int k = counts.shape[0]
    cdef int offset = 0
    cdef int i, j
    for i in range(k):
        for j in range(counts[i]):
            mutable_state[offset] = i
            offset += 1
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
def reselection(
        np.ndarray[np.float64_t, ndim=1] probs_out,
        np.ndarray[np.float64_t, ndim=1] fitnesses_in,
        np.ndarray[np.int_t, ndim=1] counts_in,
        ):
    """
    Use haploid selection.
    This corresponds to multiplicative diploid selection
    when the haplotypes are assumed to come from a diploid population.
    @param probs_out: iid child haplotype probabilities to be computed
    @param fitnesses_in: relative fitnesses of haplotypes
    @param counts_in: allele counts of the parent generation
    """
    cdef int k = counts_in.shape[0]
    cdef double x
    cdef double accum = 0
    for i in range(k):
        x = fitnesses_in[i] * counts_in[i]
        probs_out[i] = x
        accum += x
    for i in range(k):
        probs_out[i] /= accum
    return None

