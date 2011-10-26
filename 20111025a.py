"""
Eigendecompose a ue' + eu' matrix.

The vector u is a random vector
and the vector e is a constant vector of ones.
It can be proved that if u and v are both nonzero
and are not proportional to each other
then the sum of these two outer products
has exactly one positive eigenvalue and one negative eigenvalue.
"""

from StringIO import StringIO
import argparse

import numpy as np
import scipy
from scipy import linalg

import Form
import FormOut
import combobreaker

def get_form():
    """
    @return: the body of a form
    """
    # define the form objects
    form_objects = [
            Form.Integer('length', 'vector length', 5, low=2, high=10),
            Form.RadioGroup('rgroup', 'vector options', [
                Form.RadioItem(
                    'one_random',
                    'one vector is random and the other vector is constant',
                    True),
                Form.RadioItem(
                    'two_random',
                    'both vectors are random and are independent'),
                Form.RadioItem(
                    'prop',
                    'one is random and the other is a multiple of it')])]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):
    out = StringIO()
    # get a random vector and a constant vector
    u = np.random.randn(fs.length)
    if fs.one_random:
        v = np.ones(fs.length)
    elif fs.prop:
        v = u*5
    else:
        v = np.random.randn(fs.length)
    # get the symmetrized outer product
    M = np.outer(u, v) + np.outer(v, u)
    W, V = scipy.linalg.eigh(M)
    # print the results
    print >> out, 'shape %s array u:' % str(u.shape)
    print >> out, u
    print >> out
    print >> out, 'shape %s array v:' % str(v.shape)
    print >> out, v
    print >> out
    print >> out, "shape %s array M = outer(v, u) + outer(u, v):" % str(M.shape)
    print >> out, M
    print >> out
    print >> out, 'eigenvalues of M:'
    print >> out, W
    print >> out
    print >> out, 'columns are corresponding eigenvectors of M:'
    print >> out, V
    return out.getvalue()

def get_inertia(M, eps):
    inertia = {-1:0, 0:0, 1:0}
    for w in np.linalg.eigvalsh(M):
        if w < -eps:
            inertia[-1] += 1
        elif w > eps:
            inertia[1] += 1
        else:
            inertia[0] += 1
    return inertia

class TrackingChecker:
    def __init__(self, eps):
        """
        @param eps: smaller eigenvalues are said to be zero
        """
        self.eps = eps
        self.u = None
        self.v = None
        self.M = None
        self.inertia = None
    def _has_bad_inertia(self):
        return (self.inertia[-1], self.inertia[1]) != (1, 1)
    def __call__(self, uv_pair):
        self.u, self.v = uv_pair
        self.M = np.outer(self.u, self.v) + np.outer(self.v, self.u)
        self.inertia = get_inertia(self.M, self.eps)
        return self._has_bad_inertia()
    def __str__(self):
        out = StringIO()
        if self._has_bad_inertia():
            print >> out, '*** found an example with invalid inertia! ***'
            # get eigendecomposition
            W, V = scipy.linalg.eigh(self.M)
            # get response
            print >> out, 'inertia:', self.inertia
            print >> out, 'epsilon for eigenvalue zeroness:', self.eps
            print >> out, 'u:', self.u
            print >> out, 'v:', self.v
            print >> out, 'M:', self.M
            print >> out, 'eigenvalues:', W
            print >> out, 'eigenvectors as columns:'
            print >> out, V
        else:
            print >> out, 'nothing interesting was observed'
        return out.getvalue().rstrip()

def gen_uv_pairs(ndim):
    while True:
        u = np.array([1] + [0]*(ndim-1))
        v = np.random.randn(ndim)
        yield u, v

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
            '--n', type=int, default=6,
            help='use this many dimensions')
    parser.add_argument(
            '--eps', type=float, default='1e-8',
            help='smaller eigenvalues are said to be zero')
    args = parser.parse_args()
    checker = TrackingChecker(args.eps)
    print combobreaker.run_checker(checker, gen_uv_pairs(args.n))

