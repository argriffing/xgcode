"""
Are Markov processes least informative of reversible 2-state given pi and mu?

The contants pi and mu represent the stationary distribution
and something proportional to the expected number of transitions
per unit of time.
The conjecture is that Markov processes are never more informative
for divergence time than non-Markov processes with the same
summary statistics.
This web script looks for counterexamples using continuous-time hidden
Markov processes.
"""

from StringIO import StringIO
import random
import math
import itertools
from itertools import product

import numpy as np
import scipy
from scipy import linalg

import Form
import FormOut
import mrate
import ctmcmi
import msimpl
import combobreaker
import MatrixUtil
from MatrixUtil import ndot



def get_form():
    form_objects = [
            Form.Integer('nstates', 'number of states', 4, low=2, high=10),
            Form.Float('etime', 'expected divergence time',
                '2.0', low_exclusive=0, high_exclusive=10),
            ]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_mutual_information(R, A, t):
    """
    @param R: reversible Markov rate matrix
    @param A: a set of vertices on one side of the bipartition
    @param t: divergence time
    """
    n = len(R)
    A = sorted(set(A))
    B = sorted(set(range(n)) - set(A))
    v = mrate.R_to_distn(R)
    P = scipy.linalg.expm(R*t)
    J_t = (P.T * v).T
    mi = 0
    for X in (A, B):
        for Y in (A, B):
            pxy = np.sum(J_t[np.ix_(X, Y)])
            pxpy = np.sum(v[X]) * np.sum(v[Y])
            mi += pxy * math.log(pxy / pxpy)
    return mi

class Accumulate:
    def __init__(self, nstates, etime):
        self.nstates = nstates
        self.etime = etime
        self.counterexample = None
        self.n_too_close = 0
    def __call__(self):
        """
        Look for a counterexample.
        """
        # Sample a rate matrix.
        # Use a trick by Robert Kern to left and right multiply by diagonals.
        # http://mail.scipy.org/pipermail/numpy-discussion/2007-March/
        # 026809.html
        S = MatrixUtil.sample_pos_sym_matrix(self.nstates)
        v = mrate.sample_distn(self.nstates)
        R = (v**-0.5)[:,np.newaxis] * S * (v**0.5)
        R -= np.diag(np.sum(R, axis=1))
        # sample a random time
        rate = 1.0 / self.etime
        t = random.expovariate(rate)
        # sample one side of the bipartition and get the mutual information
        k = random.randrange(1, self.nstates)
        A = random.sample(range(self.nstates), k)
        mi_non_markov = get_mutual_information(R, A, t)
        # get summary statistics of the non-markov process
        Q = msimpl.get_fast_two_state(R, A)
        mi_markov = ctmcmi.get_expected_ll_ratio(Q, t)
        # check if the mutual informations are indistinguishable
        if np.allclose(mi_non_markov, mi_markov):
            self.n_too_close += 1
            return False
        if mi_non_markov < mi_markov:
            out = StringIO()
            print >> out, 'found a counterexample'
            print >> out
            print >> out, 'sampled symmetric matrix S:'
            print >> out, S
            print >> out
            print >> out, 'sampled stationary distribution v:'
            print >> out, v
            print >> out
            print >> out, 'implied rate matrix R:'
            print >> out, R
            print >> out
            print >> out, 'reduced rate matrix Q'
            print >> out, Q
            print >> out
            print >> out, 'sampled time t:', t
            print >> out
            print >> out, 'non-markov mutual information:', mi_non_markov
            print >> out, 'markov mutual information:', mi_markov
            print >> out
            self.counterexample = out.getvalue().rstrip()
            return True
    def __str__(self):
        out = StringIO()
        print >> out, 'iterations where m.i. was too close to call:',
        print >> out, self.n_too_close
        if self.counterexample:
            print >> out, self.counterexample
        else:
            print >> out, 'no counterexample was found'
        return out.getvalue().rstrip()

def get_response_content(fs):
    nseconds = 4.0
    accum = Accumulate(fs.nstates, fs.etime)
    info = combobreaker.run_callable(accum, nseconds=nseconds)
    return str(info)

