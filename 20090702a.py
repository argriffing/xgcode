"""
Check the speed of singular value decomposition.
"""

from StringIO import StringIO
import time

import numpy as np

import Form
import FormOut

g_shapes = (
        (50, 10000),
        (50, 20000),
        (50, 40000),
        (450, 10000))

def get_form():
    """
    @return: the body of a form
    """
    form_objects = []
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_matrix(n, p):
    """
    @param n: one dimension of the matrix
    @param p: another dimension of the matrix
    @return: a sampled matrix of the given dimensions
    """
    mu = 1.0
    return np.random.exponential(mu, (p, n))

def do_analysis(n, p):
    """
    @param n: one dimension of the matrix
    @param p: another dimension of the matrix
    @return: the multi-line report string
    """
    out = StringIO()
    # create the matrix
    start_time = time.time()
    M = get_matrix(n, p)
    nseconds = time.time() - start_time
    print >> out, nseconds, 'seconds',
    print >> out, 'to create a', n, 'x', p, 'matrix'
    # get the singular value decomposition
    start_time = time.time()
    np.linalg.svd(M, full_matrices=0)
    nseconds = time.time() - start_time
    print >> out, nseconds, 'seconds',
    print >> out, 'to get the singular value decomposition'
    # get the singular value decomposition of the transpose
    start_time = time.time()
    np.linalg.svd(M.T, full_matrices=0)
    nseconds = time.time() - start_time
    print >> out, nseconds, 'seconds',
    print >> out, 'to get the singular value decomposition of the transpose'
    return out.getvalue().strip()

def get_response_content(fs):
    out = StringIO()
    for n, p in g_shapes[:3]:
        print >> out, do_analysis(n, p)
        print >> out
    return out.getvalue()

def main():
    for n, p in g_shapes:
        print do_analysis(n, p)
        print

if __name__ == '__main__':
    main()

