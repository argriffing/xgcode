"""
Check the eigendecomposition of a matrix constructed from a kernel function.

The kernel function is apparently the characteristic function
of the hyperbolic cosecant distribution.
It is f(t) = t / sinh(t).
The idea is that if you feed it a vector v then
you can get a matrix Mij = (vj - vi) / sinh(vj - vi)
such that M is guaranteed to be postive definite
as long as v was not something like zero.
"""

from StringIO import StringIO
import argparse
import math

import numpy as np
import scipy
from scipy import linalg

import Form
import FormOut

def get_form():
    """
    @return: the body of a form
    """
    # define the form objects
    form_objects = [
            Form.Integer('nstates', 'number of states', 4, low=2, high=12)]
    return form_objects

def get_form_out():
    return FormOut.Report()

def sample_distribution(n):
    """
    @param n: number of states
    """
    # Get a nonnegative vector.
    v = np.random.rand(n)
    # Divide the vector by its nonnegative sum.
    distn = v / np.sum(v)
    return distn

def get_response_content(fs):
    np.set_printoptions(linewidth=200)
    out = StringIO()
    # do the analysis
    n = fs.nstates
    pi_m = sample_distribution(n)
    pi_s = sample_distribution(n)
    v = np.log(np.sqrt(pi_m / pi_s))
    K = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            x = v[j] - v[i]
            if x:
                K[i, j] = x / math.sinh(x)
            else:
                K[i, j] = 1.0
    W, V = scipy.linalg.eigh(K)
    # write the report
    print >> out, 'mutation process stationary distribution:'
    print >> out, pi_m
    print >> out
    print >> out, 'selection process stationary distribution:'
    print >> out, pi_s
    print >> out
    print >> out, 'vector to which the kernel function is applied:'
    print >> out, v
    print >> out
    print >> out, 'kernel matrix K:'
    print >> out, K
    print >> out
    print >> out, 'eigenvalues of K:'
    print >> out, W
    print >> out
    print >> out, 'eigenvectors of K:'
    print >> out, V
    return out.getvalue()

