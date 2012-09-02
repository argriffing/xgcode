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

