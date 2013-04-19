"""
Check a diagonalization of a structured three-state rate matrix.
"""

from StringIO import StringIO
import random
import math

import numpy as np

import Form
import FormOut

def get_form():
    """
    @return: the body of a form
    """
    return [
            ]

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):

    # set up print options
    np.set_printoptions(
            linewidth=1000000,
            threshold=1000000,
            )
    out = StringIO()

    # sample the random rates
    a = random.expovariate(1)
    b = random.expovariate(1)
    c = random.expovariate(1)

    # rate matrix
    Q = np.array([
        [-a, a, 0],
        [b, -(b+c), c],
        [0, 0, 0],
        ], dtype=float)

    # Compute the square root of a thing that is non-negative
    # when all rates are non-negative.
    # The non-negativity could be proved by factoring
    # the difference of squares and then using the
    # inequality of arithmetic and geometric means.
    d = math.sqrt((a+b+c)*(a+b+c) - 4*a*c)

    # According to wolfram alpha this is the diagonalization Q = U diag(w) V
    # where U V = I.
    U = np.array([
        [1, (-a+b+c-d)/(2*b), (-a+b+c+d)/(2*b)],
        [1, 1, 1],
        [1, 0, 0],
        ], dtype=float)
    w = np.array([
        0,
        -(a+b+c+d)/2,
        -(a+b+c-d)/2,
        ], dtype=float)
    V = np.array([
        [0, 0, 1],
        [-b/d, (-a+b+c+d)/(2*d), (a+b-c-d)/(2*d)],
        [b/d, (a-b-c+d)/(2*d), (-a-b+c-d)/(2*d)],
        ], dtype=float)

    # report stuff
    print >> out, 'rate a:', a
    print >> out, 'rate b:', b
    print >> out, 'rate c:', c
    print >> out
    print >> out, 'U:'
    print >> out, U
    print >> out
    print >> out, 'w:',
    print >> out, w
    print >> out
    print >> out, 'V:'
    print >> out, V
    print >> out
    print >> out, 'U V (should be the 3x3 identity):'
    print >> out, np.dot(U, V)
    print >> out
    print >> out, 'Q (the rate matrix):'
    print >> out, Q
    print >> out
    print >> out, 'U diag(w) V (should be the same as Q):'
    print >> out, np.dot(U * w, V)
    print >> out

    # show the result
    return out.getvalue()

