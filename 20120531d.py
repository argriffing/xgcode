r"""
Plot mutual information in a way that highlights the numerical difficulties.

The mutation process is site-independent 3-site 2-state-per-site.
Selection changes the relative probabilities of
even and odd parity joint states.
Look at how the mutual information at time t responds to
changes in the selection-induced log probability ratio.
"""

from StringIO import StringIO

import numpy as np

import Form
import FormOut
import ctmcmi
import RUtil
from RUtil import mk_call_str
import evozoo

def get_form():
    """
    @return: the body of a form
    """
    return [
            Form.Float('t', 'divergence time', '0.066', low_exclusive=0),
            Form.FloatInterval(
                'x_min', 'x_max', 'log probability ratio interval',
                '-2.5', '2.5'),
            Form.ImageFormat(),
            ]

def get_form_out():
    return FormOut.Image('plot')

def get_presets():
    return [
            Form.Preset(
                'expanded domain',
                {
                    't' : '0.066',
                    'x_min' : '-10',
                    'x_max' : '10',
                    'imageformat' : 'png'})]

def get_response_content(fs):
    f_info = ctmcmi.get_mutual_info_known_distn
    # define the R table headers
    headers = ['log.probability.ratio', 'mutual.information']
    # make the array
    arr = []
    for x in np.linspace(fs.x_min, fs.x_max, 101):
        row = [x]
        proc = evozoo.AlternatingHypercube_d_1(3)
        X = np.array([x])
        distn = proc.get_distn(X)
        Q = proc.get_rate_matrix(X)
        info = f_info(Q, distn, fs.t)
        row.append(info)
        arr.append(row)
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
    print >> out, mk_call_str('require', '"ggplot2"')
    print >> out, mk_call_str(
            'ggplot', 'my.table',
            mk_call_str(
                'aes',
                x='log.probability.ratio',
                y='mutual.information')),
    print >> out, '+',
    print >> out, 'geom_line()'
    return out.getvalue()

