r"""
Seek the first critical divtime for Fisher info for a 3-state process.

The mutation process is assumed to be site-independent
with three sites and two states per site.
At short times, the process without selection
gives more Fisher information.
At medium times, the process with a cycle-like lethality pattern
give more information.
"""

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
                '0.65', '0.70', low_inclusive=0),
            ]

def get_form_out():
    return FormOut.Report('results')

class Opt:
    """
    This is for brentq root finding.
    """
    def __init__(self):
        d = 3
        # precompute the hypercube
        zoo_obj = evozoo.Hypercube_d_0(d)
        self.hypercube_Q = zoo_obj.get_rate_matrix()
        self.hypercube_distn = zoo_obj.get_distn()
        # precompute the coil
        zoo_obj = evozoo.Coil_d_0(d)
        self.coil_Q = zoo_obj.get_rate_matrix()
        self.coil_distn = zoo_obj.get_distn()
    def __call__(self, t):
        """
        @param t: divergence time
        @return: signed difference between target and mutual information
        """
        hypercube_info = divtime.get_fisher_info_known_distn_fast(
                self.hypercube_Q, self.hypercube_distn, t)
        coil_info = divtime.get_fisher_info_known_distn_fast(
                self.coil_Q, self.coil_distn, t)
        return coil_info - hypercube_info

def get_response_content(fs):
    opt = Opt()
    if opt(fs.a) * opt(fs.b) >= 0:
        raise ValueError(
                'please choose a bracket that more obviously '
                'contains the right entropy')
    t_opt = scipy.optimize.brentq(opt, fs.a, fs.b)
    return str(t_opt)

