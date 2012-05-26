"""
Check formulas for information criteria of F81 processes.

For an N-state reversible Markov process,
the various information criteria are written as sums of N^2 terms,
whereas the F81 constraints lead to formulas that are sums of only N terms.
Here we numerically check that they are the same.
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
import divtime
import combobreaker
import MatrixUtil
from MatrixUtil import ndot


def get_form():
    form_objects = [
            Form.Integer('nstates', 'number of states', 4, low=2, high=10),
            ]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_f81_fi(r, v, t):
    """
    Get the Fisher information for a Felsenstein 1981 Markov process.
    @param r: randomization rate
    @param v: stationary distribution
    @param t: time
    @return: the fisher information
    """
    x = math.exp(-r*t)
    coeff = (r*x)**2 / (1 - x)
    return coeff * np.sum(v*(1-v) / (x + (1-x)*v))

def get_f81_mi(r, v, t):
    """
    Get the mutual information for a Felsenstein 1981 Markov process.
    @param r: randomization rate
    @param v: stationary distribution
    @param t: time
    @return: the mutual information
    """
    h = np.sum(v*(1-v))
    x = math.exp(-r*t)
    u = v + (1-v)*x
    return np.sum(v*u*np.log(u/v)) + h*(1-x)*math.log(1-x)

def get_f81_pollock(r, v, t):
    """
    Get the negative identity slope for a Felsenstein 1981 Markov process.
    @param r: randomization rate
    @param v: stationary distribution
    @param t: time
    @return: negative of the derivative of the expected identity proportion
    """
    h = np.sum(v*(1-v))
    x = math.exp(-r*t)
    return r*h*x

class Contradiction(Exception): pass

class Accumulate:
    def __init__(self, nstates):
        self.nstates = nstates
        self.counterexample = None
    def __call__(self):
        """
        Look for a counterexample.
        """
        n = self.nstates
        # sample a random rate and time and stationary distribution
        r = random.expovariate(1)
        t = random.expovariate(1)
        v = np.random.exponential(1, n)
        v /= np.sum(v)
        # construct the F81 rate matrix
        R = r * np.outer(np.ones(n), v)
        R -= np.diag(np.sum(R, axis=1))
        # get some information criterion values
        mi_general = ctmcmi.get_mutual_information(R, t)
        fi_general = divtime.get_fisher_information(R, t)
        mi_f81 = get_f81_mi(r, v, t)
        fi_f81 = get_f81_fi(r, v, t)
        # check for contradictions
        try:
            if not np.allclose(mi_general, mi_f81):
                raise Contradiction('mutual information')
            if not np.allclose(fi_general, fi_f81):
                raise Contradiction('fisher information')
        except Contradiction as e:
            out = StringIO()
            print >> out, 'found', str(e), 'contradiction'
            print >> out
            print >> out, 'GTR mutual information:'
            print >> out, mi_general
            print >> out
            print >> out, 'F81 mutual information:'
            print >> out, mi_f81
            print >> out
            print >> out, 'GTR Fisher information:'
            print >> out, fi_general
            print >> out
            print >> out, 'F81 Fisher information:'
            print >> out, fi_f81
            print >> out
            self.counterexample = out.getvalue()
            return True
        return False
    def __str__(self):
        out = StringIO()
        if self.counterexample:
            print >> out, self.counterexample
        else:
            print >> out, 'no counterexample was found'
        return out.getvalue().rstrip()

def get_response_content(fs):
    nseconds = 4.0
    accum = Accumulate(fs.nstates)
    info = combobreaker.run_callable(accum, nseconds=nseconds)
    return str(info)

