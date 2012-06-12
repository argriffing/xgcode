r"""
Narrow down the times after which selection can increase mutual information.

The mutation process is assumed to be site-independent
with two states per site.
If d is the number of sites
then adding selection can give mutual information near (d-1)*log(2)
regardless of divergence time.
The purpose of this web script is to determine how long it
takes the mutual information of the neutral process
to decay to this value.
"""

from StringIO import StringIO
import math

import numpy as np
import scipy
from scipy import optimize

import Form
import FormOut
import ctmcmi
import evozoo

def get_form():
    """
    @return: the body of a form
    """
    return [
            Form.Integer('d', 'number of sites', 2, low=2, high=5),
            Form.FloatInterval(
                'a', 'b', 'divtime interval',
                '.1', '.2', low_inclusive=0),
            ]

def get_form_out():
    return FormOut.Report('results')

def get_presets():
    return [
            Form.Preset(
                'three sites',
                {'d' : 3, 'a' : '0.04', 'b' : '0.10'}),
            ]

class Opt:
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

def get_response_content(fs):
    opt = Opt(fs.d)
    if opt(fs.a) * opt(fs.b) >= 0:
        raise ValueError(
                'please choose a bracket that more obviously '
                'contains the right entropy')
    t_opt = scipy.optimize.brentq(opt, fs.a, fs.b)
    return str(t_opt)

