"""
Look at properties of an effective parent independent process.

Sample a GTR rate matrix.
Then construct a constrained parent independent process rate matrix
defined by properties of the GTR rate matrix.
Check that the mutual information of the parent independent process
at time t is a lower bound of the mutual information of the
sampled process.
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
import ctmcmi
import combobreaker
import MatrixUtil
from MatrixUtil import ndot


def get_form():
    form_objects = [
            Form.Integer('nstates', 'number of states', 4, low=2, high=10)]
    return form_objects

def get_form_out():
    return FormOut.Report()

class Accumulate:
    def __init__(self, nstates):
        self.nstates = nstates
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
        # Construct a parent-independent process
        # with the same max rate and stationary distribution
        # as the sampled process.
        #max_rate = max(-np.diag(R))
        #expected_rate = np.dot(v, -np.diag(R))
        #logical_entropy = np.dot(v, 1-v)
        #randomization_rate = expected_rate / logical_entropy
        Q = np.outer(np.ones(self.nstates), v)
        Q -= np.diag(np.sum(Q, axis=1))
        Q *= max(np.diag(R) / np.diag(Q))
        # sample a random time
        t = random.expovariate(1)
        # Check that the mutual information of the
        # parent independent process is smaller.
        mi_R = ctmcmi.get_expected_ll_ratio(R, t)
        mi_Q = ctmcmi.get_expected_ll_ratio(Q, t)
        if np.allclose(mi_R, mi_Q):
            self.n_too_close += 1
            return False
        if mi_R < mi_Q:
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
            print >> out, 'parent independent process Q:'
            print >> out, Q
            print >> out
            print >> out, 'sampled time t:', t
            print >> out
            print >> out, 'mutual information of sampled process:', mi_R
            print >> out, 'mutual information of p.i. process:', mi_Q
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
    # get the user data
    nstates = fs.nstates
    # request a limited amount of time
    nseconds = 4.0
    # get the results
    accum = Accumulate(nstates)
    info = combobreaker.run_callable(accum, nseconds=nseconds)
    return str(info)

