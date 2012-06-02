r"""
Seek a critical divtime for Fisher info for a 2-state process.

The mutation process is assumed to be site-independent
with two sites and two states per site.
At medium times,
a non-lethal selection process with a single degree of freedom
gives more information.
At long times,
a path-like selection process with a single degree of freedom
gives more information.
"""

from StringIO import StringIO
import math

import numpy as np
import scipy
from scipy import optimize

import Form
import FormOut
import divtime
import evozoo

def get_form():
    """
    @return: the body of a form
    """
    return [
            Form.Float('a', 'a divtime bracket endpoint',
                '0.85', low_inclusive=0),
            Form.Float('b', 'a divtime bracket endpoint',
                '0.87', low_inclusive=0),
            ]

def get_form_out():
    return FormOut.Report('results')

class OptMinCube:
    """
    This is for minimization.
    """
    def __init__(self, t):
        self.t = t
        d = 2
        self.obj = evozoo.AlternatingHypercube_d_1(d)
    def __call__(self, X):
        Q = self.obj.get_rate_matrix(X)
        distn = self.obj.get_distn(X)
        info = divtime.get_fisher_info_known_distn_fast(Q, distn, self.t)
        return -info

class OptMinSnake:
    """
    This is for minimization.
    """
    def __init__(self, t):
        self.t = t
        d = 2
        self.obj = evozoo.AlternatingSnake_d_1(d)
    def __call__(self, X):
        Q = self.obj.get_rate_matrix(X)
        distn = self.obj.get_distn(X)
        info = divtime.get_fisher_info_known_distn_fast(Q, distn, self.t)
        return -info

class OptRoot:
    """
    This is for brentq root finding.
    """
    def __call__(self, t):
        """
        @param t: divergence time
        @return: signed difference between infos
        """
        # get the max parameterized cube info
        opt_min = OptMinCube(t)
        X0 = np.random.randn(1)
        xopt = scipy.optimize.fmin(opt_min, X0)
        max_cube_info = -opt_min(xopt)
        # get the max parameterized snake-in-the-box info
        opt_min = OptMinSnake(t)
        X0 = np.random.randn(1)
        xopt = scipy.optimize.fmin(opt_min, X0)
        max_snake_info = -opt_min(xopt)
        # return the info difference
        return max_snake_info - max_cube_info

def get_response_content(fs):
    opt = OptRoot()
    if opt(fs.a) * opt(fs.b) >= 0:
        raise ValueError(
                'please choose a bracket that more obviously '
                'contains a root')
    t_opt = scipy.optimize.brentq(opt, fs.a, fs.b)
    return str(t_opt)

