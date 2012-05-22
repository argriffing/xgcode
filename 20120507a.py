"""
Look for K-L divergence monotonicity counterexamples for CTMC processes.

K-L is Kullback-Leibler divergence,
and CTMC is Continuous Time Markov Chain.
The K-L divergence between the joint distribution of
two points at stationarity separated by time t
and the joint distribution implied by independence
(this is the mutual information between the points)
is known to decrease monotonically as the divergence time t increases.
But I am not so sure if this monotonicity is still true
when the initial state has no uncertainty;
maybe the monotonicity proof requires stationarity.
The Markov process is assumed to be irreducible and reversible.
"""

from StringIO import StringIO
import random
import math

import numpy as np
import scipy
from scipy import linalg

import Form
import FormOut
import mrate
import MatrixUtil
import combobreaker


def get_form():
    form_objects = [
            Form.Integer('nstates', 'number of states', 4, low=2, high=10),
            ]
    return form_objects

def get_form_out():
    return FormOut.Report()

class Accumulate:
    def __init__(self, nstates):
        self.nstates = nstates
        self.counterexample = None
        self.nrows_total = 0
        self.nrows_allclose = 0
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
        # Sample a short time and a longer time.
        # For each row of the transition matrix,
        # look at the K-L divergence to the stationary distribution,
        # and check that it is smaller at the larger time.
        ta, tb = sorted((random.expovariate(1), random.expovariate(1)))
        Pa = scipy.linalg.expm(R * ta)
        Pb = scipy.linalg.expm(R * tb)
        for rowa, rowb in zip(Pa, Pb):
            self.nrows_total += 1
            if np.allclose(rowa, rowb):
                self.nrows_allclose += 1
                continue
            kla = sum(x*math.log(x/y) for x, y in zip(rowa, v))
            klb = sum(x*math.log(x/y) for x, y in zip(rowb, v))
            if kla < klb:
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
                print >> out, 'sampled time ta:', ta
                print >> out, 'sampled time tb:', tb
                print >> out
                print >> out, 'transition matrix Pa:'
                print >> out, Pa
                print >> out
                print >> out, 'transition matrix Pb:'
                print >> out, Pb
                print >> out
                print >> out, 'relevant row of Pa:'
                print >> out, rowa
                print >> out
                print >> out, 'relevant row of Pb:'
                print >> out, rowb
                print >> out
                print >> out, 'K-L divergence of row of Pa from v:'
                print >> out, kla
                print >> out
                print >> out, 'K-L divergence of row of Pb from v:'
                print >> out, klb
                print >> out
                self.counterexample = out.getvalue().rstrip()
                return True
    def __str__(self):
        out = StringIO()
        print >> out, 'number of rows pairs considered:',
        print >> out, self.nrows_total
        print >> out, 'number of row pairs too similar to distinguish:',
        print >> out, self.nrows_allclose
        if self.counterexample:
            print >> out, self.counterexample
        else:
            print >> out, 'no counterexample was found'
        return out.getvalue().rstrip()

def get_response_content(fs):
    # get the user data
    nstates = fs.nstates
    # request a limited amount of time
    nseconds = 4.0
    # get the results
    accum = Accumulate(nstates)
    info = combobreaker.run_callable(accum, nseconds=nseconds)
    return str(info)

