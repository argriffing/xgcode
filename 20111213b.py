"""
Try to repeatedly improve relaxation time by adding selection.
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
        R = mrate.to_gtr_halpern_bruno(self.M, v_new)
        if not np.allclose(v_new, mrate.R_to_distn(R)):
            print v_new
            print mrate.R_to_distn(R)
            raise ValueError('stationary distribution error')
        r_sel = mrate.R_to_relaxation_time(R)
        # we want to minimize this
        return self.r_mut - r_sel

def get_response_content(fs):
    np.set_printoptions(linewidth=200)
    out = StringIO()
    n = fs.nstates
    t = 0.001
    # sample the initial mutation rate matrix
    S = sample_symmetric_rate_matrix(n)
    v = sample_distribution(n)
    M = mrate.to_gtr_halpern_bruno(S, v)
    if not np.allclose(v, mrate.R_to_distn(M)):
        raise ValueError('stationary distribution error')
    print >> out, 't:', t
    print >> out
    print >> out, 'initial GTR matrix:'
    print >> out, M
    print >> out
    # Try to iteratively increase the relaxation time
    # by repeatedly applying Halpern-Bruno selection.
    R = M
    v_old = v
    for i in range(20):
        # print some properties of the matrix
        print >> out, v_old
        print >> out, mrate.R_to_relaxation_time(R)
        print >> out
        f = MyOpt(R, t)
        x0 = [1.0] * (n - 1)
        result = scipy.optimize.fmin(
                f, x0, disp=0, full_output=1, ftol=0.000001)
        xopt, fopt, niters, funcalls, warnflag = result
        if fopt > 0:
            print >> out, 'failed to increase relaxation time'
            print >> out
            break
        # compute the next stationary distribution
        v_target = X_to_distn(xopt)
        v_new = (1 - t) * v_old + t * v_target
        print >> out, v_new - v_old
        print >> out
        # compute the next rate matrix and update its stationary distribution
        R = mrate.to_gtr_halpern_bruno(R, v_new)
        if not np.allclose(v_new, mrate.R_to_distn(R)):
            raise ValueError('stationary distribution error')
        v_old = v_new
    print >> out, 'final rate matrix:'
    print >> out, R
    print >> out
    return out.getvalue()

