r"""
Look for general selection that maximizes mutual information.

The mutation process is assumed to be a
site-independent d-site 2-state-per-site process
with uniform stationary distribution.
The selection is assumed to work according to a formula used
by Sella and Hirsh and Halpern and Bruno.
No other constraint is applied to the selection.
"""

from StringIO import StringIO
import argparse
import math
import time
import random
from itertools import product

import numpy as np
import scipy
from scipy import linalg
from scipy import optimize

import Form
import FormOut
import ctmcmi
import mrate
import divtime
import cheeger
import MatrixUtil
from MatrixUtil import ndot
import evozoo


def get_form():
    """
    @return: the body of a form
    """
    return [
            Form.Integer('d', 'number of sites', 2, low=2, high=5),
            Form.Float('a', 'a divtime bracket endpoint',
                '.1', low_inclusive=0),
            Form.Float('b', 'a divtime bracket endpoint',
                '.2', low_inclusive=0),
            ]

def get_form_out():
    return FormOut.Report('results')

def get_presets():
    return [
            Form.Preset(
                'three sites',
                {'d' : 3, 'a' : '0.04', 'b' : '0.10'}),
            ]

class OptRoot:
    """
    This is for brentq root finding.
    """
    def __init__(self, d):
        """
        @param d: number of sites
        """
        self.d = d
        self.target = (d-1)*math.log(2)
        # precompute the matrices
        zoo_obj = evozoo.Hypercube_d_0(d)
        self.distn = zoo_obj.get_distn()
        self.Q = zoo_obj.get_rate_matrix()
    def __call__(self, t):
        """
        @param t: divergence time
        @return: signed difference between target and mutual information
        """
        info_value = ctmcmi.get_mutual_info_known_distn_fast(
                self.Q, self.distn, t)
        return self.target - info_value

class OptMin:
    """
    This is for minimization.
    """
    def __init__(self, d, t):
        """
        @param d: number of sites
        @param t: divergence time
        """
        self.t = t
        self.zoo_obj = evozoo.GeneralHypercube_d_full(d)
    def __call__(self, X):
        """
        @param X: log probability ratios
        @return: negative of mutual information
        """
        Q = self.zoo_obj.get_rate_matrix(X)
        distn = self.zoo_obj.get_distn(X)
        return -ctmcmi.get_mutual_info_known_distn_fast(Q, distn, self.t)

def get_opt_divtime(d, a, b):
    opt_root = OptRoot(d)
    if opt_root(a) * opt_root(b) >= 0:
        raise ValueError(
                'please choose a bracket that more obviously '
                'contains the right entropy')
    return scipy.optimize.brentq(opt_root, a, b)

def do_search_opt(opt_min, nseconds):
    """
    @param opt_min: a function object for numerical optimization
    @param nseconds: spend roughly this many seconds searching
    @return: a freeform multiline string report
    """
    np.set_printoptions(linewidth=200)
    # init the search
    df = opt_min.zoo_obj.get_df()
    best_info = None
    best_xopt = None
    t0 = time.time()
    niter = 0
    while time.time() - t0 < nseconds:
        niter += 1
        X0 = np.random.randn(df)
        # fmin_powell is not so great and also slow
        # fmin_bfgs seems pretty good and also fast
        # well sometimes fmin_bfgs is taking forever
        # but maybe this is only with the high variance energies
        xopt = scipy.optimize.fmin(
                opt_min, X0, maxiter=10000, maxfun=10000)
        """
        xopt = scipy.optimize.fmin_bfgs(opt_min, X0)
        """
        """
        xopt = scipy.optimize.brute(
                opt_min, ranges=tuple([(-10, 10, 5)]*df))
        """
        info = -opt_min(xopt)
        if best_info is None or info > best_info:
            best_info = info
            best_xopt = xopt
    # write the report
    out = StringIO()
    Q = opt_min.zoo_obj.get_rate_matrix(best_xopt)
    distn = opt_min.zoo_obj.get_distn(best_xopt)
    print >> out, 'numerical maximization runs:', niter
    print >> out, 't:', opt_min.t
    print >> out, 'expected info:', (opt_min.zoo_obj.d - 1)*math.log(2)
    print >> out, 'observed info:', best_info
    print >> out, distn
    print >> out, Q
    return out.getvalue().rstrip()

def get_response_content(fs):
    nseconds = 4
    t = get_opt_divtime(fs.d, fs.a, fs.b)
    print t
    opt_min = OptMin(fs.d, t)
    return do_search_opt(opt_min, nseconds) + '\n'

