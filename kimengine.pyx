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
def kimura_integrand(float x, float c, float d):
    cdef float n2cx = -2.*c*x
    cdef float retvalue =  exp(n2cx*d*(1.-x) + n2cx)
    return retvalue

