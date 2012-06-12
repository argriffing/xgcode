r"""
Seek the second critical divtime for Fisher info for a 3-state process.

The mutation process is assumed to be site-independent
with three sites and two states per site.
At medium times, a process with a cycle-like lethality pattern
gives more information.
At long times, a non-lethal selection process
parameterized by a single degree of freedom
gives more information.
"""

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
            Form.FloatInterval(
                'a', 'b', 'divtime interval',
                '1.13', '1.15', low_inclusive=0),
            ]

def get_form_out():
    return FormOut.Report('results')

class OptMin:
    """
    This is for minimization.
    """
    def __init__(self, t):
        self.t = t
        d = 3
        self.obj = evozoo.AlternatingHypercube_d_1(d)
    def __call__(self, X):
        Q = self.obj.get_rate_matrix(X)
        distn = self.obj.get_distn(X)
        info = divtime.get_fisher_info_known_distn_fast(Q, distn, self.t)
        return -info

class OptRoot:
    """
    This is for brentq root finding.
    """
    def __init__(self):
        d = 3
        # precompute the coil
        zoo_obj = evozoo.Coil_d_0(d)
        self.coil_Q = zoo_obj.get_rate_matrix()
        self.coil_distn = zoo_obj.get_distn()
    def __call__(self, t):
        """
        @param t: divergence time
        @return: signed difference between infos
        """
        # get the max parameterized cube info
        opt_min = OptMin(t)
        X0 = np.random.randn(1)
        xopt = scipy.optimize.fmin(opt_min, X0)
        max_cube_info = -opt_min(xopt)
        # get the coil info
        coil_info = divtime.get_fisher_info_known_distn_fast(
                self.coil_Q, self.coil_distn, t)
        # return the info difference
        return max_cube_info - coil_info

def get_response_content(fs):
    opt = OptRoot()
    if opt(fs.a) * opt(fs.b) >= 0:
        raise ValueError(
                'please choose a bracket that more obviously '
                'contains a root')
    t_opt = scipy.optimize.brentq(opt, fs.a, fs.b)
    return str(t_opt)

