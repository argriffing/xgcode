"""
Fast functions for Kimura integrals.

This is related to population genetics.
For compilation instructions see
http://docs.cython.org/src/reference/compilation.html
For example:
$ cython -a kimengine.pyx
$ gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing \
      -I/usr/include/python2.7 -o kimengine.so kimengine.c
"""

import numpy as np
cimport numpy as np
cimport cython
from libc.math cimport log, exp

np.import_array()

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def kimura_integrand(double x, double c, double d):
    cdef double n2cx = -2.*c*x
    cdef double retvalue =  exp(n2cx*d*(1.-x) + n2cx)
    return retvalue

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def denom_poly(double c, double d):
    """
    This is a polynomial approximation for small c and small d.
    But the smallness of c is more important.
    """
    cdef double accum = 0.
    cdef int k = 0
    #
    k = 0
    cdef double a_0 = 1.
    cdef double b_0 = (1./1.) * 1.
    accum += a_0 * b_0
    #
    k = 1
    cdef double a_1 = -a_0 * (2.*c) / k
    cdef double b_1 = (1./6.) * (3. + d)
    accum += a_1 * b_1
    #
    k = 2
    cdef double a_2 = -a_1 * (2.*c) / k
    cdef double b_2 = (1./30.) * (10. + d*(5. + d))
    accum += a_2 * b_2
    #
    k = 3
    cdef double a_3 = -a_2 * (2.*c) / k
    cdef double b_3 = (1./140.) * (35. + d*(21. + d*(7. + d)))
    accum += a_3 * b_3
    #
    k = 4
    cdef double a_4 = -a_3 * (2.*c) / k
    cdef double b_4 = (1./630.) * (126. + d*(84. + d*(36. + d*(9. + d))))
    accum += a_4 * b_4
    #
    k = 5
    cdef double a_5 = -a_4 * (2.*c) / k
    cdef double b_5 = (1./2772.) * (
            462. + d*(330. + d*(165. + d*(55. + d*(11. + d)))))
    accum += a_5 * b_5
    #
    k = 6
    cdef double a_6 = -a_5 * (2.*c) / k
    cdef double b_6 = (1./12012.) * (
            1716. + d*(1287. + d*(715. + d*(286. + d*(78. + d*(13. + d))))))
    accum += a_6 * b_6
    #
    return accum

