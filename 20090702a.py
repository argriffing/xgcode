"""Check the speed of singular value decomposition.
"""

from StringIO import StringIO
import random
import math
import time

import numpy

from SnippetUtil import HandlingError
import MatrixUtil
import Form
import NewickIO

def get_form():
    """
    @return: the body of a form
    """
    form_objects = []
    return form_objects

def get_matrix(n, p):
    """
    @param n: one dimension of the matrix
    @param p: another dimension of the matrix
    @return: a sampled matrix of the given dimensions
    """
    M = numpy.zeros((p, n))
    for i in range(p):
        for j in range(n):
            M[i,j] = random.expovariate(1.0)
    return M

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
    print >> out, nseconds, 'seconds to create a', n, 'x', p, 'matrix'
    # get the singular value decomposition
    start_time = time.time()
    numpy.linalg.svd(M, full_matrices=0)
    nseconds = time.time() - start_time
    print >> out, nseconds, 'seconds to get the singular value decomposition'
    # get the singular value decomposition of the transpose
    start_time = time.time()
    numpy.linalg.svd(M.T, full_matrices=0)
    nseconds = time.time() - start_time
    print >> out, nseconds, 'seconds to get the singular value decomposition of the transpose'
    return out.getvalue().strip()

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    out = StringIO()
    # get the response
    n = 50
    for p in (10000, 20000, 40000):
        report = do_analysis(n, p)
        print >> out, report
        print >> out
    n = 30*15
    report = do_analysis(n, 10000)
    print >> out, report
    print >> out

    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
