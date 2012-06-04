r"""
Plot max Fisher information for several selection models.

The mutation process is site-independent d-site 2-state-per-site.
The site-dependent selection processes vary, but they all fall under
the Halpern-Bruno-like formula that associates a set of selection parameters
with each state.
Because this is a numerical optimization in possibly high dimensional space,
the search will often fail to find a global max information.
"""

from StringIO import StringIO

import numpy as np
import scipy
from scipy import optimize

import Form
import FormOut
import divtime
import evozoo
import cbreaker
import progrid
import RUtil
from RUtil import mk_call_str


# variable name, description, python object
g_process_triples = [
        ('general_cube', 'hypercube with general selection',
            evozoo.GeneralHypercube_d_full),
        ('one_param_cube', 'hypercube with single parameter stationary distn',
            evozoo.AlternatingHypercube_d_1),
        ('distinguished_corner_pair', 'a corner pair is distinguished',
            evozoo.DistinguishedCornerPairHypercube_d_1),
        ('one_param_cycle', 'max induced cycle with 1-parameter distn',
            evozoo.AlternatingCoil_d_1),
        ('one_param_path', 'max induced path with 1-parameter distn',
            evozoo.AlternatingSnake_d_1),
        ('cube', 'hypercube with uniform stationary distribution',
            evozoo.Hypercube_d_0),
        ('cycle', 'maximal induced cycle',
            evozoo.Coil_d_0),
        ('path', 'maximal induced path',
            evozoo.Snake_d_0),
        ]

def get_form():
    """
    @return: the body of a form
    """
    check_items = [Form.CheckItem(a, b, True) for a, b, c in g_process_triples]
    return [
            Form.Integer('d', 'number of sites', 3, low=2, high=5),
            Form.CheckGroup(
                'processes',
                'plot max divtime info for these parameterized processes',
                check_items),
            Form.Float('start_time', 'start time', '0.6', low_exclusive=0),
            Form.Float('stop_time', 'stop time', '1.4', low_exclusive=0),
            Form.ImageFormat(),
            ]

def get_form_out():
    return FormOut.Image('plot')

class OptDep:
    def __init__(self, zoo_obj, t, f_info):
        """
        @param zoo_obj: an object from the evozoo module
        @param t: divergence time
        @param f_info: info function that takes a rate matrix and a time
        """
        self.zoo_obj = zoo_obj
        self.t = t
        self.f_info = f_info
    def __call__(self, X):
        """
        @param X: some log ratio probabilities
        @return: neg info value for minimization
        """
        distn = self.zoo_obj.get_distn(X)
        Q = self.zoo_obj.get_rate_matrix(X)
        return -self.f_info(Q, distn, self.t)

def get_response_content(fs):
    # validate and store user input
    if fs.stop_time <= fs.start_time:
        raise ValueError('check the start and stop times')
    f_info = divtime.get_fisher_info_known_distn_fast
    requested_triples = []
    for triple in g_process_triples:
        name, desc, zoo_obj = triple
        if getattr(fs, name):
            requested_triples.append(triple)
    if not requested_triples:
        raise ValueError('nothing to plot')
    # define the R table headers
    r_names = [a.replace('_', '.') for a, b, c in requested_triples]
    headers = ['t'] + r_names
    # Spend a lot of time doing the optimizations
    # to construct the points for the R table.
    arr = []
    for t in cbreaker.throttled(
            progrid.gen_binary(fs.start_time, fs.stop_time),
            nseconds=5, ncount=200):
        row = [t]
        for python_name, desc, zoo_class in requested_triples:
            zoo_obj = zoo_class(fs.d)
            df = zoo_obj.get_df()
            opt_dep = OptDep(zoo_obj, t, f_info)
            if df:
                X0 = np.random.randn(df)
                xopt = scipy.optimize.fmin(
                        opt_dep, X0, maxiter=10000, maxfun=10000)
                # I would like to use scipy.optimize.minimize
                # except that this requires a newer version of
                # scipy than is packaged for ubuntu right now.
                # fmin_bfgs seems to have problems sometimes
                # either hanging or maxiter=10K is too big.
                """
                xopt = scipy.optimize.fmin_bfgs(opt_dep, X0,
                        gtol=1e-8, maxiter=10000)
                """
            else:
                xopt = np.array([])
            info_value = -opt_dep(xopt)
            row.append(info_value)
        arr.append(row)
    arr.sort()
    npoints = len(arr)
    # create the R table string and scripts
    # get the R table
    table_string = RUtil.get_table_string(arr, headers)
    # get the R script
    script = get_ggplot()
    # create the R plot image
    device_name = Form.g_imageformat_to_r_function[fs.imageformat]
    retcode, r_out, r_err, image_data = RUtil.run_plotter(
            table_string, script, device_name)
    if retcode:
        raise RUtil.RError(r_err)
    return image_data

def get_ggplot():
    out = StringIO()
    print >> out, mk_call_str('require', '"reshape"')
    print >> out, mk_call_str('require', '"ggplot2"')
    print >> out, 'my.table.long <-',
    print >> out, mk_call_str('melt', 'my.table', id='"t"')
    print >> out, 'ggplot(data=my.table.long,'
    print >> out, mk_call_str('aes', x='t', y='value', colour='variable')
    print >> out, ') + geom_line()',
    return out.getvalue()

