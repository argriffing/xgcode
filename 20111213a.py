"""
Check the effect of cleverly chosen selection on rate matrix relaxation time.
"""

from StringIO import StringIO
import argparse
import math

import numpy as np
import scipy
from scipy import linalg

import Form
import FormOut
from MatrixUtil import ndot
import mrate
import divtime
import combobreaker


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

def sample_symmetric_rate_matrix(n):
    """
    @param n: number of states
    """
    # Get a nonnegative asymmetric matrix.
    M = np.random.rand(n, n)
    # Symmetrize by adding to the transpose.
    S = M + M.T
    # Subtract row sum from diagonal.
    R = S - np.diag(np.sum(S, axis=1))
    return R

def X_to_distn(X):
    """
    @param X: a vector to be converted into a finite distribution
    @return: a finite distribution as a numpy array
    """
    v_target = np.exp(list(X) + [1])
    return v_target / np.sum(v_target)

class MyOpt:
    def __init__(self, M, t):
        """
        @param M: mutation matrix
        @param t: the distance to go in the requested direction
        """
        self.M = M
        self.t = t
        # get the stationary distribution of the mutation process
        self.v = mrate.R_to_distn(M)
        # get the mutation process relaxation time
        self.r_mut = mrate.R_to_relaxation_time(M)
    def __call__(self, X):
        """
        @param X: a vector to be converted into a finite distribution
        """
        v_target = X_to_distn(X)
        v_new = (1 - self.t) * self.v + self.t * v_target
        R = divtime.to_gtr_halpern_bruno(self.M, v_new)
        if not np.allclose(v_new, mrate.R_to_distn(R)):
            raise ValueError('stationary distribution error')
        r_sel = mrate.R_to_relaxation_time(R)
        # we want to minimize this
        return self.r_mut - r_sel

def maxind(arr):
    value, index = max((v, i) for i, v in enumerate(arr))
    return index

class Checker:
    def __init__(self, nstates):
        self.nstates = nstates
        #self.t = 0.001
        self.t = 0.001
        self.M = None
        self.fiedler = None
        self.r_mut = None
        self.r_sel = None
    def _get_opt_target(self):
        f = MyOpt(self.M, self.t)
        x0 = [1.0] * (self.nstates - 1)
        result = scipy.optimize.fmin(
                f, x0, disp=0, full_output=1, ftol=0.000001)
        xopt, fopt, niters, funcalls, warnflag = result
        print fopt
        if fopt < 0:
            return X_to_distn(xopt)
        else:
            return None
    def __call__(self):
        """
        @return: True if a counterexample is found
        """
        n = self.nstates
        # sample a fairly generic GTR mutation rate matrix
        S = sample_symmetric_rate_matrix(n)
        v = sample_distribution(n)
        M = divtime.to_gtr_halpern_bruno(S, v)
        # look at the fiedler-like eigenvector of the mutation rate matrix
        r_recip, fiedler = mrate._R_to_eigenpair(M)
        r_mut = 1 / r_recip
        value_min, state_min = min((fiedler[i], i) for i in range(n))
        value_max, state_max = max((fiedler[i], i) for i in range(n))
        # move the stationary distribution towards a 50/50 distribution
        v_target = np.zeros(n)
        v_target[state_min] = 0.5
        v_target[state_max] = 0.5
        v_new = (1 - self.t) * v + self.t * v_target
        R = divtime.to_gtr_halpern_bruno(M, v_new)
        r_sel = mrate.R_to_relaxation_time(R)
        # the mutation-selection balance should have longer relaxation time
        #if r_sel < r_mut:
        #if True:
        if maxind(np.abs(fiedler / v)) != maxind(np.abs(fiedler / np.sqrt(v))):
            self.M = M
            self.fiedler = fiedler
            self.r_mut = r_mut
            self.r_sel = r_sel
            self.v = v
            self.v_new = v_new
            self.v_target = v_target
            self.opt_target = self._get_opt_target()
            return True
        else:
            return False
    def __str__(self):
        out = StringIO()
        print >> out, 'mutation rate matrix:'
        print >> out, self.M
        print >> out
        print >> out, 'mutation process fiedler-like vector:'
        print >> out, self.fiedler
        print >> out
        #print >> out, 'mutation process relaxation time:'
        #print >> out, self.r_mut
        #print >> out
        #print >> out, 'mutation-selection process relaxation time:'
        #print >> out, self.r_sel
        #print >> out
        print >> out, 'mutation process stationary distribution:'
        print >> out, self.v
        print >> out
        print >> out, 'abs of fiedler divided by stationary probability:'
        print >> out, np.abs(self.fiedler / self.v)
        print >> out
        print >> out, 'abs of fiedler divided by sqrt stationary probability:'
        print >> out, np.abs(self.fiedler / np.sqrt(self.v))
        print >> out
        print >> out, 'numeric-solver-based target:'
        print >> out, self.opt_target
        print >> out
        print >> out, 'fiedler-based target:'
        print >> out, self.v_target
        print >> out
        print >> out, 'mutation-selection process stationary distribution:'
        print >> out, self.v_new
        print >> out
        return out.getvalue()


def process(nstates, nseconds):
    f = Checker(nstates)
    info = combobreaker.run_callable(f, nseconds=nseconds, niterations=None)
    return str(info)

def get_response_content(fs):
    np.set_printoptions(linewidth=200)
    nseconds = 5.0
    return process(fs.nstates, nseconds)
